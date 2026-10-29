# %% imports
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
PROJECT_DIR = project_dir("phase4")
STAND_DATA_FILEPATH_PIRKANMAA = (
    inputs_dir_for_project(PROJECT_DIR) / "stand_data_pirkanmaa.json"
)
STAND_DATA_FILEPATH_POHJANMAA = (
    inputs_dir_for_project(PROJECT_DIR) / "stand_data_pohjois_pohjanmaa.json"
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
    - Spruce is the dominant species of a site with sfc 3
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

    # Removing spruces as dominant species in sites with fertility class = 2 or 3
    if stand_data.allometry_file_per_layer[
        "dominant"
    ].species_id == 2 and stand_data.site_fertility_class in {2, 3}:
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


def predict_ditch_depth(drainage_age, peat_thickness=0.61, ditch_bed_slope=0.62):
    """
    Predict ditch shallowing by Hökkä et al. 2020, Baltic Forestry 26(2), article id 453. https://doi.org/10.46490/BF453
    """
    base = (
        49.139
        + 23.981 * np.power(drainage_age, -0.2)
        - 0.343 * drainage_age
        - 0.404 * drainage_age * peat_thickness
        + 17.961 * peat_thickness
        + 3.224 * ditch_bed_slope
    )
    if drainage_age > 35:
        base += 0.00134 * np.power(drainage_age, 2)
    return base


def estimate_drainage_age(ditch_depth_cm, peat_thickness=0.61, ditch_bed_slope=0.62):
    # Ditch depth difference at a given age
    def ditch_depth_difference(age):
        return (
            predict_ditch_depth(age, peat_thickness, ditch_bed_slope) - ditch_depth_cm
        )

    # Try solving in a reasonable range of drainage ages
    solution = root_scalar(ditch_depth_difference, bracket=[1, 100], method="brentq")

    if solution.converged:
        return np.round(solution.root)
    else:
        return None


def get_ditch_shallowing(ditch_depth, time=20):
    if ditch_depth <= -0.27:
        est_drainage_age = estimate_drainage_age((-1) * 100 * ditch_depth)
        ditch_shallowing = predict_ditch_depth(est_drainage_age + time)
        return np.round((-1) * 0.01 * ditch_shallowing, 2)
    else:
        return round(ditch_depth * 0.8, 2)


def create_simulation_parameters(
    maakunta: str,
    stand_data_document: StandDataDocument,
) -> list[SimulationParams]:

    simulation_parameters: list[SimulationParams] = []

    N_SOIL_COLS = 20
    # TODO: fix years
    INITIAL_YEAR = 2004
    END_YEAR = 2020

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
            BASE_SCENARIOS = ("control", "DNM_80", "ash_fertilization")

            for scenario_name in BASE_SCENARIOS:
                ditch_depth = -0.60
                if scenario_name == "DNM_85":
                    ditch_depth = -0.85

                fertilization_params = None
                if scenario_name == "ash_fertilization":
                    fertilization_params = AshFertilizationParameters(
                        application_year=INITIAL_YEAR + 10,
                        fertilizer_dose=5000,
                        K_in_ash=100,
                        P_in_ash=45,
                    )

                simulation_parameters.append(
                    SimulationParams(
                        metadata=SimulationMetaData(
                            project_dir=PROJECT_DIR,
                            run_id="run_01",
                            stand_id=stand_id,
                            scenario_id=scenario_name,
                            tags={"region": maakunta},
                        ),
                        susi_params=SusiParams(
                            weather_parameters=WeatherParams(
                                FMI_weather_filepath=data_dir_for_project(PROJECT_DIR)
                                / "weather"
                                / "CFw.csv",
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
                                location=LocationsForPhotoParams(
                                    "Hyytiala"
                                )  # Change for Oulu area!
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
                                ditch_depth_west=[ditch_depth],
                                ditch_depth_east=[ditch_depth],
                                ditch_depth_20y_west=[
                                    get_ditch_shallowing(ditch_depth, time=20)
                                ],
                                ditch_depth_20y_east=[
                                    get_ditch_shallowing(ditch_depth, time=20)
                                ],
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
                                fertilization=fertilization_params,
                                peat_temperature=PeatTemperatureParams(),
                            ),
                        ),
                    )
                )
    return simulation_parameters


# %%

# Finally, create the list of parameters that will go into the susi simulation
pirkanmaa_parameters = create_simulation_parameters(
    maakunta="pirkanmaa", stand_data_document=stand_data_document_pirkanmaa
)
pohjanmaa_parameters = create_simulation_parameters(
    maakunta="pohjanmaa", stand_data_document=stand_data_document_pohjanmaa
)

# TODO: fix this
all_parameters = pirkanmaa_parameters[:30] + pohjanmaa_parameters[:30]


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
