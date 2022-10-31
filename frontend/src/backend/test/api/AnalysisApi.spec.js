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

(function (root, factory) {
  if (typeof define === "function" && define.amd) {
    // AMD.
    define(["expect.js", process.cwd() + "/src/index"], factory);
  } else if (typeof module === "object" && module.exports) {
    // CommonJS-like environments that support module.exports, like Node.
    factory(require("expect.js"), require(process.cwd() + "/src/index"));
  } else {
    // Browser globals (root is window)
    factory(root.expect, root.Ccfatigue);
  }
})(this, function (expect, Ccfatigue) {
  "use strict";

  var instance;

  beforeEach(function () {
    instance = new Ccfatigue.AnalysisApi();
  });

  var getProperty = function (object, getter, property) {
    // Use getter method if present; otherwise, get the property directly.
    if (typeof object[getter] === "function") return object[getter]();
    else return object[property];
  };

  var setProperty = function (object, setter, property, value) {
    // Use setter method if present; otherwise, set the property directly.
    if (typeof object[setter] === "function") object[setter](value);
    else object[property] = value;
  };

  describe("AnalysisApi", function () {
    describe("runCldFileAnalysisCldFilePost", function () {
      it("should call runCldFileAnalysisCldFilePost successfully", function (done) {
        //uncomment below and update the code to test runCldFileAnalysisCldFilePost
        //instance.runCldFileAnalysisCldFilePost(function(error) {
        //  if (error) throw error;
        //expect().to.be();
        //});
        done();
      });
    });
    describe("runCycleCountingFileAnalysisCycleCountingFilePost", function () {
      it("should call runCycleCountingFileAnalysisCycleCountingFilePost successfully", function (done) {
        //uncomment below and update the code to test runCycleCountingFileAnalysisCycleCountingFilePost
        //instance.runCycleCountingFileAnalysisCycleCountingFilePost(function(error) {
        //  if (error) throw error;
        //expect().to.be();
        //});
        done();
      });
    });
    describe("runDamageSummationFileAnalysisDamageSummationFilePost", function () {
      it("should call runDamageSummationFileAnalysisDamageSummationFilePost successfully", function (done) {
        //uncomment below and update the code to test runDamageSummationFileAnalysisDamageSummationFilePost
        //instance.runDamageSummationFileAnalysisDamageSummationFilePost(function(error) {
        //  if (error) throw error;
        //expect().to.be();
        //});
        done();
      });
    });
    describe("runFatigueFailureFileAnalysisFatigueFailureFilePost", function () {
      it("should call runFatigueFailureFileAnalysisFatigueFailureFilePost successfully", function (done) {
        //uncomment below and update the code to test runFatigueFailureFileAnalysisFatigueFailureFilePost
        //instance.runFatigueFailureFileAnalysisFatigueFailureFilePost(function(error) {
        //  if (error) throw error;
        //expect().to.be();
        //});
        done();
      });
    });
    describe("runSnCurveFileAnalysisSnCurveFilePost", function () {
      it("should call runSnCurveFileAnalysisSnCurveFilePost successfully", function (done) {
        //uncomment below and update the code to test runSnCurveFileAnalysisSnCurveFilePost
        //instance.runSnCurveFileAnalysisSnCurveFilePost(function(error) {
        //  if (error) throw error;
        //expect().to.be();
        //});
        done();
      });
    });
  });
});
