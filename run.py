# %% imports
from typing import Sequence
from dataclasses import dataclass, replace
from susi.io.load_output_data import StandID
import datetime
import argparse
from multiprocessing import Pool
import numpy as np
from scipy.optimize import root_scalar

from susi.core.susi_main import Susi
from susi.io.execution_config import MultipleSusis, SimulationParams
from susi.io.metadata_model import SimulationMetaData

from susi.io.project_layout import (
    data_dir_for_project,
    project_dir,
    inputs_dir_for_project,
)
from susi.io.susi_parameter_model import (
    CanopyLayerAllometry,
    CanopyLayerName,
    CanopyParams,
    LocationsForPhotoParams,
    OrganicLayerParams,
    OutputParams,
    PeatTemperatureParams,
    PeatTypes,
    SimulationConfig,
    SiteParams,
    StandParams,
    SusiParams,
    WeatherParams,
    get_photo_parameters_by_location,
    h_mor_from_drainage_and_mass_mor_Pitkanen,
    AshFertilizationParameters,
)
from susi.io.stand_data import (
    load_stand_data_document_from_json,
    StandDataDocument,
)


# %% Parse CLI arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    "n_parallel_processes", help="Number of parallel processes", type=int
)
cli_args = parser.parse_args()

# %% Constants
PROJECT_DIR = project_dir("wq_rosa")
STAND_DATA_FILEPATH_PIRKANMAA = (
    inputs_dir_for_project(PROJECT_DIR) / "stand_data_pirkanmaa.json"
)
STAND_DATA_FILEPATH_POHJANMAA = (
    inputs_dir_for_project(PROJECT_DIR) / "stand_data_pohjanmaa.json"
)

# %% Load stand data
stand_data_document_pirkanmaa = load_stand_data_document_from_json(
    STAND_DATA_FILEPATH_PIRKANMAA
)
stand_data_document_pohjanmaa = load_stand_data_document_from_json(
    STAND_DATA_FILEPATH_POHJANMAA
)

# %% Inspect and filter stand data

stands_pirkanmaa = stand_data_document_pirkanmaa.stands

MAP_SPECIES_NUMBER = {1: "pine", 2: "spruce", 3: "other"}
SITE_FERTILITY_CLASSES = (2, 3, 4, 5)


def stands_per_sfc(stands):
    sfc_counts = {i: 0 for i in SITE_FERTILITY_CLASSES}
    species_per_sfc_counts = {
        i: {"pine": 0, "spruce": 0, "other": 0} for i in SITE_FERTILITY_CLASSES
    }
    for stand_id, stand_data in stands.items():
        sfc = stand_data.site_fertility_class
        species_number_dominant = stand_data.allometry_file_per_layer[
            "dominant"
        ].species_id
        species_name_dominant = MAP_SPECIES_NUMBER[species_number_dominant]
        sfc_counts[sfc] += 1
        species_per_sfc_counts[sfc][species_name_dominant] += 1

    print("Site fertility class counts:")
    print(species_per_sfc_counts)


print("Pirkanmaa")
stands_per_sfc(stand_data_document_pirkanmaa.stands)

print("\nP.Pohjanmaa")
stands_per_sfc(stand_data_document_pohjanmaa.stands)

print("Conclusion:")
print("Do not use data for sfc =2, as there are too few stands")
print("Do not use data where dominant species=3")


def build_stand_params(
    stand_data_document: StandDataDocument,
    stand_id: StandID,
    n_cols: int,
) -> StandParams | None:
    """
    Builds the StandParams for a single stand.
    Returns None if:
    - sfc=2, because there are few datapoints
    - dominant species=3 (not pine nor spruce)
    - Spruce is the dominant species of a site with sfc 4 or 5
    """
    stand_data = stand_data_document.stands[stand_id]

    canopy_layer_allometry = CanopyLayerAllometry.with_single_allometry_per_layer(
        layers=stand_data.allometry_file_per_layer,
        n=n_cols,
    )
    if stand_data.site_fertility_class == 2:
        return None

    # Only pine and spruce can be the dominant species
    if stand_data.allometry_file_per_layer["dominant"].species_id == 3:
        return None

    # Removing spruces as dominant species in sites with fertility class = 4 or 5
    if stand_data.allometry_file_per_layer[
        "dominant"
    ].species_id == 2 and stand_data.site_fertility_class in {4, 5}:
        return None

    return StandParams(
        site_fertility_class=stand_data.site_fertility_class,
        canopy_layer_allometry=canopy_layer_allometry,
        # The engine needs an age for every layer, including ones the stand
        # has no allometry for: those start at 0.0.
        initial_canopylayer_age_years={
            layer: stand_data.initial_age_per_layer.get(layer, 0.0)
            for layer in CanopyLayerName
        },
    )


# %% Create the scenarios

INITIAL_YEAR = 1975
END_YEAR = 2025


@dataclass(frozen=True)
class Scenario:
    name: str
    ditch_depth_east: float
    ditch_depth_west: float
    ditch_depth_20y_east: float
    ditch_depth_20y_west: float
    fertilization: AshFertilizationParameters | None


def create_scenarios() -> list[Scenario]:
    base_scenarios: list[Scenario] = [
        Scenario(
            name="DD60-60_NS",
            ditch_depth_east=-0.6,
            ditch_depth_west=-0.6,
            ditch_depth_20y_east=-0.6,
            ditch_depth_20y_west=-0.6,
            fertilization=None,
        ),
        Scenario(
            name="DD60-60_S",
            ditch_depth_east=-0.6,
            ditch_depth_west=-0.6,
            ditch_depth_20y_east=-0.3,
            ditch_depth_20y_west=-0.3,
            fertilization=None,
        ),
        Scenario(
            name="DD90-90_S",
            ditch_depth_east=-0.9,
            ditch_depth_west=-0.9,
            ditch_depth_20y_east=-0.6,
            ditch_depth_20y_west=-0.6,
            fertilization=None,
        ),
        Scenario(
            name="WM30-30_S",
            ditch_depth_east=-0.3,
            ditch_depth_west=-0.3,
            ditch_depth_20y_east=-0.2,
            ditch_depth_20y_west=-0.2,
            fertilization=None,
        ),
        Scenario(
            name="WM60-10_S",
            ditch_depth_east=-0.6,
            ditch_depth_west=-0.1,
            ditch_depth_20y_east=-0.4,
            ditch_depth_20y_west=-0.1,
            fertilization=None,
        ),
        Scenario(
            name="WM30-10_S",
            ditch_depth_east=-0.3,
            ditch_depth_west=-0.1,
            ditch_depth_20y_east=-0.2,
            ditch_depth_20y_west=-0.1,
            fertilization=None,
        ),
    ]
    # Add ash fertilization to each scenario:
    # Each base scenario gets a fertilized twin
    fertilization_params = AshFertilizationParameters(
        application_year=INITIAL_YEAR,
        fertilizer_dose=4000,
        K_in_ash=120,
        P_in_ash=45,
    )

    ash_scenarios = [
        replace(
            scenario,
            name=f"{scenario.name}_ash",
            fertilization=fertilization_params,
        )
        for scenario in base_scenarios
    ]

    return base_scenarios + ash_scenarios


def create_simulation_parameters(
    maakunta: str, stand_data_document: StandDataDocument, scenarios: list[Scenario]
) -> list[SimulationParams]:

    simulation_parameters: list[SimulationParams] = []

    N_SOIL_COLS = 20
    weather_files = {"pirkanmaa": "pirkanmaa.csv", "pohjanmaa": "pohjanmaa.csv"}

    for stand_id, stand_data in stand_data_document.stands.items():
        stand_params = build_stand_params(
            stand_data_document=stand_data_document,
            stand_id=stand_id,
            n_cols=N_SOIL_COLS,
        )
        # Remove stand if sfc=2 and/or dominant species is not pine/spruce
        if stand_params is None:
            continue
        else:
            for scenario in scenarios:
                simulation_parameters.append(
                    SimulationParams(
                        metadata=SimulationMetaData(
                            project_dir=PROJECT_DIR,
                            run_id="run_01",
                            stand_id=stand_id,
                            scenario_id=scenario.name,
                            tags={"region": maakunta},
                        ),
                        susi_params=SusiParams(
                            weather_parameters=WeatherParams(
                                FMI_weather_filepath=data_dir_for_project(PROJECT_DIR)
                                / "weather"
                                / weather_files[maakunta],
                            ),
                            simulation_config=SimulationConfig(
                                start_date=datetime.datetime(INITIAL_YEAR, 1, 1),
                                end_date=datetime.datetime(END_YEAR, 12, 31),
                            ),
                            stand_params=stand_params,
                            canopy_parameters=CanopyParams(),
                            organic_layer_parameters=OrganicLayerParams(),
                            output_parameters=OutputParams(),
                            photo_parameters=get_photo_parameters_by_location(
                                location=LocationsForPhotoParams("All_data")
                            ),
                            site_parameters=SiteParams(
                                L=40.0,
                                n=N_SOIL_COLS,
                                sitename="susirun",
                                sfc_specification=1,
                                hdom=None,
                                vol=None,
                                smc="Peatland",
                                nLyrs=60,
                                dzLyr=0.05,
                                ditch_depth_west=[scenario.ditch_depth_west],
                                ditch_depth_east=[scenario.ditch_depth_east],
                                ditch_depth_20y_west=[scenario.ditch_depth_20y_west],
                                ditch_depth_20y_east=[scenario.ditch_depth_20y_east],
                                scenario_name=["deprecated"],  # kasvunlisaykset
                                drain_age=50.0,
                                initial_h=-0.2,
                                slope=0.0,
                                peat_type=[PeatTypes.generic] * 8,
                                peat_type_bottom=[PeatTypes.generic],
                                anisotropy=10.0,
                                vonP=True,
                                vonP_top=[2, 5, 5, 5, 6, 6, 7, 7],
                                vonP_bottom=8,
                                bd_top=None,
                                bd_bottom=0.16,
                                peatN=None,
                                peatP=None,
                                peatK=None,
                                enable_peattop=True,
                                enable_peatmiddle=True,
                                enable_peatbottom=True,
                                rho_mor=90.0,
                                h_mor=h_mor_from_drainage_and_mass_mor_Pitkanen,
                                cutting_management=None,
                                depoN=4.0,
                                depoP=0.1,
                                depoK=1.0,
                                fertilization=scenario.fertilization,
                                peat_temperature=PeatTemperatureParams(),
                            ),
                        ),
                    )
                )
    return simulation_parameters


# %%

# Finally, create the list of parameters that will go into the susi simulation
scenarios = create_scenarios()
pirkanmaa_parameters = create_simulation_parameters(
    maakunta="pirkanmaa",
    stand_data_document=stand_data_document_pirkanmaa,
    scenarios=scenarios,
)
pohjanmaa_parameters = create_simulation_parameters(
    maakunta="pohjanmaa",
    stand_data_document=stand_data_document_pohjanmaa,
    scenarios=scenarios,
)

all_parameters = pirkanmaa_parameters + pohjanmaa_parameters


execution_config = MultipleSusis(
    simulation_parameter_list=all_parameters,
    n_parallel_processes=cli_args.n_parallel_processes,
)


def run_susi(simulation_parameters: SimulationParams) -> None:
    # Initiate susi class
    susi = Susi(simulation_parameters)

    # Run simulation
    susi.run()

    # Save stuff
    susi.write_params_and_metadata()


# %% Execute parallel processing
if __name__ == "__main__":
    with Pool(processes=execution_config.n_parallel_processes) as pool:
        pool.map(func=run_susi, iterable=execution_config.simulation_parameter_list)
