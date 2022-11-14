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

import ApiClient from "../ApiClient";

/**
 * The QuasiStaticTest model module.
 * @module model/QuasiStaticTest
 * @version 0.1.0
 */
class QuasiStaticTest {
  /**
   * Constructs a new <code>QuasiStaticTest</code>.
   * @alias module:model/QuasiStaticTest
   * @param crackDisplacement {Array.<Number>}
   * @param crackLoad {Array.<Number>}
   * @param crackLength {Array.<Number>}
   * @param displacement {Object.<String, Array.<Number>>}
   * @param load {Object.<String, Array.<Number>>}
   * @param strain {Object.<String, Array.<Number>>}
   * @param stress {Object.<String, Array.<Number>>}
   */
  constructor(
    crackDisplacement,
    crackLoad,
    crackLength,
    displacement,
    load,
    strain,
    stress
  ) {
    QuasiStaticTest.initialize(
      this,
      crackDisplacement,
      crackLoad,
      crackLength,
      displacement,
      load,
      strain,
      stress
    );
  }

  /**
   * Initializes the fields of this object.
   * This method is used by the constructors of any subclasses, in order to implement multiple inheritance (mix-ins).
   * Only for internal use.
   */
  static initialize(
    obj,
    crackDisplacement,
    crackLoad,
    crackLength,
    displacement,
    load,
    strain,
    stress
  ) {
    obj["crack_displacement"] = crackDisplacement;
    obj["crack_load"] = crackLoad;
    obj["crack_length"] = crackLength;
    obj["displacement"] = displacement;
    obj["load"] = load;
    obj["strain"] = strain;
    obj["stress"] = stress;
  }

  /**
   * Constructs a <code>QuasiStaticTest</code> from a plain JavaScript object, optionally creating a new instance.
   * Copies all relevant properties from <code>data</code> to <code>obj</code> if supplied or a new instance if not.
   * @param {Object} data The plain JavaScript object bearing properties of interest.
   * @param {module:model/QuasiStaticTest} obj Optional instance to populate.
   * @return {module:model/QuasiStaticTest} The populated <code>QuasiStaticTest</code> instance.
   */
  static constructFromObject(data, obj) {
    if (data) {
      obj = obj || new QuasiStaticTest();

      if (data.hasOwnProperty("crack_displacement")) {
        obj["crack_displacement"] = ApiClient.convertToType(
          data["crack_displacement"],
          ["Number"]
        );
      }
      if (data.hasOwnProperty("crack_load")) {
        obj["crack_load"] = ApiClient.convertToType(data["crack_load"], [
          "Number",
        ]);
      }
      if (data.hasOwnProperty("crack_length")) {
        obj["crack_length"] = ApiClient.convertToType(data["crack_length"], [
          "Number",
        ]);
      }
      if (data.hasOwnProperty("displacement")) {
        obj["displacement"] = ApiClient.convertToType(data["displacement"], {
          String: ["Number"],
        });
      }
      if (data.hasOwnProperty("load")) {
        obj["load"] = ApiClient.convertToType(data["load"], {
          String: ["Number"],
        });
      }
      if (data.hasOwnProperty("strain")) {
        obj["strain"] = ApiClient.convertToType(data["strain"], {
          String: ["Number"],
        });
      }
      if (data.hasOwnProperty("stress")) {
        obj["stress"] = ApiClient.convertToType(data["stress"], {
          String: ["Number"],
        });
      }
    }
    return obj;
  }
}

/**
 * @member {Array.<Number>} crack_displacement
 */
QuasiStaticTest.prototype["crack_displacement"] = undefined;

/**
 * @member {Array.<Number>} crack_load
 */
QuasiStaticTest.prototype["crack_load"] = undefined;

/**
 * @member {Array.<Number>} crack_length
 */
QuasiStaticTest.prototype["crack_length"] = undefined;

/**
 * @member {Object.<String, Array.<Number>>} displacement
 */
QuasiStaticTest.prototype["displacement"] = undefined;

/**
 * @member {Object.<String, Array.<Number>>} load
 */
QuasiStaticTest.prototype["load"] = undefined;

/**
 * @member {Object.<String, Array.<Number>>} strain
 */
QuasiStaticTest.prototype["strain"] = undefined;

/**
 * @member {Object.<String, Array.<Number>>} stress
 */
QuasiStaticTest.prototype["stress"] = undefined;

export default QuasiStaticTest;
