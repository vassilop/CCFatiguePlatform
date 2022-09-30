/**
 * ccfatigue
 * No description provided (generated by Openapi Generator https://github.com/openapitools/openapi-generator)
 *
 * The version of the OpenAPI document: 0.1.0
 *
 *
 * NOTE: This class is auto generated by OpenAPI Generator (https://openapi-generator.tech).
 * https://openapi-generator.tech
 * Do not edit the class manually.
 *
 */
/**
 * Enum class ExperimentFieldNames.
 * @enum {}
 * @readonly
 */
export default class ExperimentFieldNames {
  /**
   * value: "fracture_mode"
   * @const
   */
  fracture_mode = "fracture_mode";

  /**
   * value: "material_type_fiber_material"
   * @const
   */
  material_type_fiber_material = "material_type_fiber_material";

  /**
   * value: "material_type_resin"
   * @const
   */
  material_type_resin = "material_type_resin";

  /**
   * value: "laminates_and_assemblies_stacking_sequence"
   * @const
   */
  laminates_and_assemblies_stacking_sequence =
    "laminates_and_assemblies_stacking_sequence";

  /**
   * Returns a <code>ExperimentFieldNames</code> enum value from a Javascript object name.
   * @param {Object} data The plain JavaScript object containing the name of the enum value.
   * @return {module:model/ExperimentFieldNames} The enum <code>ExperimentFieldNames</code> value.
   */
  static constructFromObject(object) {
    return object;
  }
}
