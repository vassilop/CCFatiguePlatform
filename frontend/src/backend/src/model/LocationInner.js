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
 * The LocationInner model module.
 * @module model/LocationInner
 * @version 0.1.0
 */
class LocationInner {
  /**
   * Constructs a new <code>LocationInner</code>.
   * @alias module:model/LocationInner
   */
  constructor() {
    LocationInner.initialize(this);
  }

  /**
   * Initializes the fields of this object.
   * This method is used by the constructors of any subclasses, in order to implement multiple inheritance (mix-ins).
   * Only for internal use.
   */
  static initialize(obj) {}

  /**
   * Constructs a <code>LocationInner</code> from a plain JavaScript object, optionally creating a new instance.
   * Copies all relevant properties from <code>data</code> to <code>obj</code> if supplied or a new instance if not.
   * @param {Object} data The plain JavaScript object bearing properties of interest.
   * @param {module:model/LocationInner} obj Optional instance to populate.
   * @return {module:model/LocationInner} The populated <code>LocationInner</code> instance.
   */
  static constructFromObject(data, obj) {
    if (data) {
      obj = obj || new LocationInner();
    }
    return obj;
  }
}

export default LocationInner;
