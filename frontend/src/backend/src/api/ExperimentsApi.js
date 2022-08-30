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
import DashboardPlots from "../model/DashboardPlots";
import ExperimentDataPreprocessed from "../model/ExperimentDataPreprocessed";
import ExperimentFieldNames from "../model/ExperimentFieldNames";
import HTTPValidationError from "../model/HTTPValidationError";
import PageExperimentModel from "../model/PageExperimentModel";

/**
 * Experiments service.
 * @module api/ExperimentsApi
 * @version 0.1.0
 */
export default class ExperimentsApi {
  /**
   * Constructs a new ExperimentsApi.
   * @alias module:api/ExperimentsApi
   * @class
   * @param {module:ApiClient} [apiClient] Optional API client implementation to use,
   * default to {@link module:ApiClient#instance} if unspecified.
   */
  constructor(apiClient) {
    this.apiClient = apiClient || ApiClient.instance;
  }

  /**
   * Callback function to receive the result of the getExperimentsExperimentsGet operation.
   * @callback module:api/ExperimentsApi~getExperimentsExperimentsGetCallback
   * @param {String} error Error message, if any.
   * @param {module:model/PageExperimentModel} data The data returned by the service call.
   * @param {String} response The complete HTTP response.
   */

  /**
   * Get Experiments
   * Get all experiments
   * @param {Object} opts Optional parameters
   * @param {String} opts.query  (default to '')
   * @param {String} opts.textSearch  (default to '')
   * @param {Number} opts.page  (default to 1)
   * @param {Number} opts.size  (default to 50)
   * @param {module:api/ExperimentsApi~getExperimentsExperimentsGetCallback} callback The callback function, accepting three arguments: error, data, response
   * data is of type: {@link module:model/PageExperimentModel}
   */
  getExperimentsExperimentsGet(opts, callback) {
    opts = opts || {};
    let postBody = null;

    let pathParams = {};
    let queryParams = {
      query: opts["query"],
      text_search: opts["textSearch"],
      page: opts["page"],
      size: opts["size"],
    };
    let headerParams = {};
    let formParams = {};

    let authNames = [];
    let contentTypes = [];
    let accepts = ["application/json"];
    let returnType = PageExperimentModel;
    return this.apiClient.callApi(
      "/experiments",
      "GET",
      pathParams,
      queryParams,
      headerParams,
      formParams,
      postBody,
      authNames,
      contentTypes,
      accepts,
      returnType,
      null,
      callback
    );
  }

  /**
   * Callback function to receive the result of the getFieldDistinctExperimentsFieldDistinctGet operation.
   * @callback module:api/ExperimentsApi~getFieldDistinctExperimentsFieldDistinctGetCallback
   * @param {String} error Error message, if any.
   * @param {Array.<String>} data The data returned by the service call.
   * @param {String} response The complete HTTP response.
   */

  /**
   * Get Field Distinct
   * Get all distinct values for field column, sorted
   * @param {module:model/ExperimentFieldNames} field
   * @param {module:api/ExperimentsApi~getFieldDistinctExperimentsFieldDistinctGetCallback} callback The callback function, accepting three arguments: error, data, response
   * data is of type: {@link Array.<String>}
   */
  getFieldDistinctExperimentsFieldDistinctGet(field, callback) {
    let postBody = null;
    // verify the required parameter 'field' is set
    if (field === undefined || field === null) {
      throw new Error(
        "Missing the required parameter 'field' when calling getFieldDistinctExperimentsFieldDistinctGet"
      );
    }

    let pathParams = {
      field: field,
    };
    let queryParams = {};
    let headerParams = {};
    let formParams = {};

    let authNames = [];
    let contentTypes = [];
    let accepts = ["application/json"];
    let returnType = ["String"];
    return this.apiClient.callApi(
      "/experiments/{field}/distinct",
      "GET",
      pathParams,
      queryParams,
      headerParams,
      formParams,
      postBody,
      authNames,
      contentTypes,
      accepts,
      returnType,
      null,
      callback
    );
  }

  /**
   * Callback function to receive the result of the getTestsDashboardPlotsExperimentsTestsDashboardPlotsGet operation.
   * @callback module:api/ExperimentsApi~getTestsDashboardPlotsExperimentsTestsDashboardPlotsGetCallback
   * @param {String} error Error message, if any.
   * @param {module:model/DashboardPlots} data The data returned by the service call.
   * @param {String} response The complete HTTP response.
   */

  /**
   * Get Tests Dashboard Plots
   * Return the 4 Bokeh plots used in Test Dashboard  Note: as we don't have real data yet, we hard code things this so it will render the 10 first tests of the experiment 1 (only experiment we have) : + experiment=1 + 1<tests_ids<10 then we mascarade test_id field so that it looks like to be matching the one asked for.
   * @param {Object} opts Optional parameters
   * @param {Number} opts.experimentId
   * @param {Array.<Number>} opts.testIds
   * @param {module:api/ExperimentsApi~getTestsDashboardPlotsExperimentsTestsDashboardPlotsGetCallback} callback The callback function, accepting three arguments: error, data, response
   * data is of type: {@link module:model/DashboardPlots}
   */
  getTestsDashboardPlotsExperimentsTestsDashboardPlotsGet(opts, callback) {
    opts = opts || {};
    let postBody = null;

    let pathParams = {};
    let queryParams = {
      experiment_id: opts["experimentId"],
      test_ids: this.apiClient.buildCollectionParam(opts["testIds"], "multi"),
    };
    let headerParams = {};
    let formParams = {};

    let authNames = [];
    let contentTypes = [];
    let accepts = ["application/json"];
    let returnType = DashboardPlots;
    return this.apiClient.callApi(
      "/experiments/tests_dashboard_plots",
      "GET",
      pathParams,
      queryParams,
      headerParams,
      formParams,
      postBody,
      authNames,
      contentTypes,
      accepts,
      returnType,
      null,
      callback
    );
  }

  /**
   * Callback function to receive the result of the postDataPreprocessCheckExperimentsDataPreprocessCheckPost operation.
   * @callback module:api/ExperimentsApi~postDataPreprocessCheckExperimentsDataPreprocessCheckPostCallback
   * @param {String} error Error message, if any.
   * @param {module:model/ExperimentDataPreprocessed} data The data returned by the service call.
   * @param {String} response The complete HTTP response.
   */

  /**
   * Post Data Preprocess Check
   * @param {File} file
   * @param {module:api/ExperimentsApi~postDataPreprocessCheckExperimentsDataPreprocessCheckPostCallback} callback The callback function, accepting three arguments: error, data, response
   * data is of type: {@link module:model/ExperimentDataPreprocessed}
   */
  postDataPreprocessCheckExperimentsDataPreprocessCheckPost(file, callback) {
    let postBody = null;
    // verify the required parameter 'file' is set
    if (file === undefined || file === null) {
      throw new Error(
        "Missing the required parameter 'file' when calling postDataPreprocessCheckExperimentsDataPreprocessCheckPost"
      );
    }

    let pathParams = {};
    let queryParams = {};
    let headerParams = {};
    let formParams = {
      file: file,
    };

    let authNames = [];
    let contentTypes = ["multipart/form-data"];
    let accepts = ["application/json"];
    let returnType = ExperimentDataPreprocessed;
    return this.apiClient.callApi(
      "/experiments/data_preprocess_check",
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
      null,
      callback
    );
  }
}
