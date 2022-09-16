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

import ApiClient from "./ApiClient";
import AppInfo from "./model/AppInfo";
import DashboardPlots from "./model/DashboardPlots";
import EchartLine from "./model/EchartLine";
import ExperimentDataPreprocessed from "./model/ExperimentDataPreprocessed";
import ExperimentFieldNames from "./model/ExperimentFieldNames";
import ExperimentModel from "./model/ExperimentModel";
import HTTPValidationError from "./model/HTTPValidationError";
import LocationInner from "./model/LocationInner";
import PageExperimentModel from "./model/PageExperimentModel";
import PageTestModel from "./model/PageTestModel";
import Plots from "./model/Plots";
import SnCurveMethod from "./model/SnCurveMethod";
import SnCurveResult from "./model/SnCurveResult";
import TestModel from "./model/TestModel";
import TestPlot from "./model/TestPlot";
import UnitInfo from "./model/UnitInfo";
import ValidationError from "./model/ValidationError";
import AnalysisApi from "./api/AnalysisApi";
import DefaultApi from "./api/DefaultApi";
import ExperimentsApi from "./api/ExperimentsApi";
import TestsApi from "./api/TestsApi";

/**
 * JS API client generated by OpenAPI Generator.<br>
 * The <code>index</code> module provides access to constructors for all the classes which comprise the public API.
 * <p>
 * An AMD (recommended!) or CommonJS application will generally do something equivalent to the following:
 * <pre>
 * var Ccfatigue = require('index'); // See note below*.
 * var xxxSvc = new Ccfatigue.XxxApi(); // Allocate the API class we're going to use.
 * var yyyModel = new Ccfatigue.Yyy(); // Construct a model instance.
 * yyyModel.someProperty = 'someValue';
 * ...
 * var zzz = xxxSvc.doSomething(yyyModel); // Invoke the service.
 * ...
 * </pre>
 * <em>*NOTE: For a top-level AMD script, use require(['index'], function(){...})
 * and put the application logic within the callback function.</em>
 * </p>
 * <p>
 * A non-AMD browser application (discouraged) might do something like this:
 * <pre>
 * var xxxSvc = new Ccfatigue.XxxApi(); // Allocate the API class we're going to use.
 * var yyy = new Ccfatigue.Yyy(); // Construct a model instance.
 * yyyModel.someProperty = 'someValue';
 * ...
 * var zzz = xxxSvc.doSomething(yyyModel); // Invoke the service.
 * ...
 * </pre>
 * </p>
 * @module index
 * @version 0.1.0
 */
export {
  /**
   * The ApiClient constructor.
   * @property {module:ApiClient}
   */
  ApiClient,
  /**
   * The AppInfo model constructor.
   * @property {module:model/AppInfo}
   */
  AppInfo,
  /**
   * The DashboardPlots model constructor.
   * @property {module:model/DashboardPlots}
   */
  DashboardPlots,
  /**
   * The EchartLine model constructor.
   * @property {module:model/EchartLine}
   */
  EchartLine,
  /**
   * The ExperimentDataPreprocessed model constructor.
   * @property {module:model/ExperimentDataPreprocessed}
   */
  ExperimentDataPreprocessed,
  /**
   * The ExperimentFieldNames model constructor.
   * @property {module:model/ExperimentFieldNames}
   */
  ExperimentFieldNames,
  /**
   * The ExperimentModel model constructor.
   * @property {module:model/ExperimentModel}
   */
  ExperimentModel,
  /**
   * The HTTPValidationError model constructor.
   * @property {module:model/HTTPValidationError}
   */
  HTTPValidationError,
  /**
   * The LocationInner model constructor.
   * @property {module:model/LocationInner}
   */
  LocationInner,
  /**
   * The PageExperimentModel model constructor.
   * @property {module:model/PageExperimentModel}
   */
  PageExperimentModel,
  /**
   * The PageTestModel model constructor.
   * @property {module:model/PageTestModel}
   */
  PageTestModel,
  /**
   * The Plots model constructor.
   * @property {module:model/Plots}
   */
  Plots,
  /**
   * The SnCurveMethod model constructor.
   * @property {module:model/SnCurveMethod}
   */
  SnCurveMethod,
  /**
   * The SnCurveResult model constructor.
   * @property {module:model/SnCurveResult}
   */
  SnCurveResult,
  /**
   * The TestModel model constructor.
   * @property {module:model/TestModel}
   */
  TestModel,
  /**
   * The TestPlot model constructor.
   * @property {module:model/TestPlot}
   */
  TestPlot,
  /**
   * The UnitInfo model constructor.
   * @property {module:model/UnitInfo}
   */
  UnitInfo,
  /**
   * The ValidationError model constructor.
   * @property {module:model/ValidationError}
   */
  ValidationError,
  /**
   * The AnalysisApi service constructor.
   * @property {module:api/AnalysisApi}
   */
  AnalysisApi,
  /**
   * The DefaultApi service constructor.
   * @property {module:api/DefaultApi}
   */
  DefaultApi,
  /**
   * The ExperimentsApi service constructor.
   * @property {module:api/ExperimentsApi}
   */
  ExperimentsApi,
  /**
   * The TestsApi service constructor.
   * @property {module:api/TestsApi}
   */
  TestsApi,
};
