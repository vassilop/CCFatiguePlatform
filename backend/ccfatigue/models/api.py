"""
API Models
"""

from enum import Enum
from pydantic import BaseModel
from typing import Optional


class OrmModel(BaseModel):
    """
    Parent class for all API models
    """

    class Config:
        orm_mode = True


class ExperimentModel(OrmModel):
    """
    Defines how experiment is seen on the API
    """

    id: int
    laboratory: Optional[str]
    researcher: str
    date: Optional[str]
    experiment_type: str

    fracture: bool
    fracture_mode: Optional[str]

    fatigue_test_type: Optional[str]
    quasi_static_test_type: Optional[str]
    temperature_test_type: Optional[str]

    measuring_equipment: Optional[str]
    reliability_level: Optional[float]

    control_mode: Optional[str]

    publication_title: Optional[str]
    publication_author: Optional[str]
    publication_year: Optional[str]
    publication_doi: Optional[str]
    publication_images_repository: Optional[str]

    material_type_sample_type: Optional[str]
    material_type_fiber_material: Optional[str]
    material_type_fiber_form: Optional[str]
    material_type_area_density: Optional[float]
    material_type_resin: Optional[str]
    material_type_hardener: Optional[str]
    material_type_mixing_ratio: Optional[str]

    laminates_and_assemblies_curing_time: Optional[float]
    laminates_and_assemblies_curing_temperature: Optional[float]
    laminates_and_assemblies_curing_pressure: Optional[float]
    laminates_and_assemblies_fiber_volume_ratio: Optional[float]
    laminates_and_assemblies_stacking_sequence: Optional[str]

    measurement_measuring_points: Optional[int]

    dic_analysis_subset_size: Optional[int]
    dic_analysis_step_size: Optional[int]


class TestModel(OrmModel):
    """
    Defines how test is seen on the API
    """

    id: int

    experiment_id: int

    specimen_number: Optional[str]
    stress_ratio: Optional[float]
    maximum_stress: Optional[float]
    loading_rate: Optional[float]
    run_out: Optional[bool]


class ExperimentFieldNames(str, Enum):
    """
    Enumerate all fields that can be requested as distinct
    """

    fracture_mode = "fracture_mode"
    material_type_fiber_material = "material_type_fiber_material"
    material_type_resin = "material_type_resin"
    laminates_and_assemblies_stacking_sequence = (
        "laminates_and_assemblies_stacking_sequence"
    )
