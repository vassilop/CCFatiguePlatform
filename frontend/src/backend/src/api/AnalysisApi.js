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
import SnCurveResult from "../model/SnCurveResult";

/**
 * Analysis service.
 * @module api/AnalysisApi
 * @version 0.1.0
 */
export default class AnalysisApi {
  /**
   * Constructs a new AnalysisApi.
   * @alias module:api/AnalysisApi
   * @class
   * @param {module:ApiClient} [apiClient] Optional API client implementation to use,
   * default to {@link module:ApiClient#instance} if unspecified.
   */
  constructor(apiClient) {
    this.apiClient = apiClient || ApiClient.instance;
  }

  /**
   * Run Cycle Counting File
   * @param {module:model/CycleCountingMethod} method
   * @param {File} file
   * @return {Promise} a {@link https://www.promisejs.org/|Promise}, with an object containing data of type {@link File} and HTTP response
   */
  runCycleCountingFileAnalysisCycleCountingFilePostWithHttpInfo(method, file) {
    let postBody = null;
    // verify the required parameter 'method' is set
    if (method === undefined || method === null) {
      throw new Error(
        "Missing the required parameter 'method' when calling runCycleCountingFileAnalysisCycleCountingFilePost"
      );
    }
    // verify the required parameter 'file' is set
    if (file === undefined || file === null) {
      throw new Error(
        "Missing the required parameter 'file' when calling runCycleCountingFileAnalysisCycleCountingFilePost"
      );
    }

    let pathParams = {};
    let queryParams = {
      method: method,
    };
    let headerParams = {};
    let formParams = {
      file: file,
    };

    let authNames = [];
    let contentTypes = ["multipart/form-data"];
    let accepts = ["application/json"];
    let returnType = File;
    return this.apiClient.callApi(
      "/analysis/cycleCounting/file",
      "POST",
      pathParams,
      queryParams,
      headerParams,
      formParams,
      postBody,
      authNames,
      contentTypes,
      accepts,
      returnType,
      null
    );
  }

  /**
   * Run Cycle Counting File
   * @param {module:model/CycleCountingMethod} method
   * @param {File} file
   * @return {Promise} a {@link https://www.promisejs.org/|Promise}, with data of type {@link File}
   */
  runCycleCountingFileAnalysisCycleCountingFilePost(method, file) {
    return this.runCycleCountingFileAnalysisCycleCountingFilePostWithHttpInfo(
      method,
      file
    ).then(function (response_and_data) {
      return response_and_data.data;
    });
  }

  /**
   * Run Sn Curve File
   * @param {Array.<module:model/SnCurveMethod>} methods
   * @param {Array.<Number>} rRatios
   * @param {File} file
   * @return {Promise} a {@link https://www.promisejs.org/|Promise}, with an object containing data of type {@link module:model/SnCurveResult} and HTTP response
   */
  runSnCurveFileAnalysisSnCurveFilePostWithHttpInfo(methods, rRatios, file) {
    let postBody = null;
    // verify the required parameter 'methods' is set
    if (methods === undefined || methods === null) {
      throw new Error(
        "Missing the required parameter 'methods' when calling runSnCurveFileAnalysisSnCurveFilePost"
      );
    }
    // verify the required parameter 'rRatios' is set
    if (rRatios === undefined || rRatios === null) {
      throw new Error(
        "Missing the required parameter 'rRatios' when calling runSnCurveFileAnalysisSnCurveFilePost"
      );
    }
    // verify the required parameter 'file' is set
    if (file === undefined || file === null) {
      throw new Error(
        "Missing the required parameter 'file' when calling runSnCurveFileAnalysisSnCurveFilePost"
      );
    }

    let pathParams = {};
    let queryParams = {
      methods: this.apiClient.buildCollectionParam(methods, "multi"),
      rRatios: this.apiClient.buildCollectionParam(rRatios, "multi"),
    };
    let headerParams = {};
    let formParams = {
      file: file,
    };

    let authNames = [];
    let contentTypes = ["multipart/form-data"];
    let accepts = ["application/json"];
    let returnType = SnCurveResult;
    return this.apiClient.callApi(
      "/analysis/snCurve/file",
      "POST",
      pathParams,
      queryParams,
      headerParams,
      formParams,
      postBody,
      authNames,
      contentTypes,
      accepts,
      returnType,
      null
    );
  }

  /**
   * Run Sn Curve File
   * @param {Array.<module:model/SnCurveMethod>} methods
   * @param {Array.<Number>} rRatios
   * @param {File} file
   * @return {Promise} a {@link https://www.promisejs.org/|Promise}, with data of type {@link module:model/SnCurveResult}
   */
  runSnCurveFileAnalysisSnCurveFilePost(methods, rRatios, file) {
    return this.runSnCurveFileAnalysisSnCurveFilePostWithHttpInfo(
      methods,
      rRatios,
      file
    ).then(function (response_and_data) {
      return response_and_data.data;
    });
  }
}
