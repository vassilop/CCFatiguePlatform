export default {
  namespaced: true,
  state: {
    filteredExperiments: {
      filters: {}, // filters in UI
      experiments: [], // the experiments fetched from API
      pagination: {
        page: 1, // page number of current fetched experiments. Overwrirten by UI
        size: 0, // max number of fetched experiments returned per request. Overwrirten by UI
        total: 0, // number of experiments in DB matching the filters. Given by backend
      },
      loading: false, // used to let the user know that data is loading
    },
    allFractureMode: [],
    allMaterialTypeFiberMaterial: [],
    allMaterialTypeResin: [],
    allLaminatesAndAssembliesStackingSequence: [],
    oneExperiment: {
      experimentId: null, // id of experiment
      experiment: {}, // the experiment itself
      tests: [], // list of tests fetched from API for that experiment
      pagination: {
        page: 1, // page number of current fetched tests. Overwrirten by UI
        size: 0, // max number of fetched tests returned per request. Overwrirten by UI
        total: 0, // number of tests in DB matching the filters. Given by backend
      },
      loading: false, // used to let the user know that data is loading
      loadingTests: false,
    },
    units: {},
  },
  mutations: {
    nowWeLoadFilteredExperiments(state, payload) {
      state.filteredExperiments = {
        filters: payload.filters,
        experiments: [],
        pagination: {
          page: payload.pagination.page,
          size: payload.pagination.size,
          total: 0,
        },
        loading: true,
      };
    },
    nowWeLoadOneExperiment(state, payload) {
      state.oneExperiment = {
        experimentId: payload.experimentId,
        experiment: {},
        tests: [],
        pagination: {
          page: payload.pagination.page,
          size: payload.pagination.size,
          total: 0,
        },
        loading: true,
        loadingTests: true,
      };
    },
    storeFilteredExperiments(state, data) {
      state.filteredExperiments = {
        ...state.filteredExperiments,
        experiments: data.items,
        pagination: {
          page: data.page,
          size: data.size,
          total: data.total,
        },
        loading: false,
      };
    },
    emptyFilteredExperiments(state) {
      state.filteredExperiments = {
        ...state.filteredExperiments,
        experiments: [],
        pagination: {
          page: 1,
          size: state.filteredExperiments.pagination.size,
          total: 0,
        },
        loading: false,
      };
    },
    storeAllFractureMode(state, data) {
      state.allFractureMode = data;
    },
    storeAllMaterialTypeFiberMaterial(state, data) {
      state.allMaterialTypeFiberMaterial = data;
    },
    storeAllMaterialTypeResin(state, data) {
      state.allMaterialTypeResin = data;
    },
    storeAllLaminatesAndAssembliesStackingSequence(state, data) {
      state.allLaminatesAndAssembliesStackingSequence = data;
    },
    storeOneExperiment(state, data) {
      state.oneExperiment = {
        ...state.oneExperiment,
        experiment: data.items[0],
        loading: false,
      };
    },
    storeOneExperimentTests(state, data) {
      state.oneExperiment = {
        ...state.oneExperiment,
        tests: data.items,
        pagination: {
          page: data.page,
          size: data.size,
          total: data.total,
        },
        loadingTests: false,
      };
    },
    storeUnits(state, data) {
      state.units = data.reduce((acc, item) => {
        acc[item.subject] = item.unit;
        return acc;
      }, {});
    },
  },
  actions: {
    fetchUnits({ commit, state }) {
      if (Object.keys(state.units).length === 0) {
        this._vm.$defaultApi.getUnitsUnitsGet((error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeUnits", data);
          }
        });
      }
    },
    fetchOneExperimentWithTests({ commit, state }, payload) {
      commit("nowWeLoadOneExperiment", payload);
      this._vm.$experimentsApi.getExperimentsExperimentsGet(
        { query: `id:${payload.experimentId}` },
        (error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeOneExperiment", data);
          }
        }
      );
      this._vm.$testsApi.getTestsTestsGet(
        payload.experimentId,
        {
          page: state.oneExperiment.pagination.page,
          size: state.oneExperiment.pagination.size,
        },
        (error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeOneExperimentTests", data);
          }
        }
      );
    },
    fetchFilteredExperiments({ commit, state }, payload) {
      commit("nowWeLoadFilteredExperiments", payload);

      const queryElements = [];

      // type (FA|QS)
      const types = [
        payload.filters.typeFA ? "FA" : null,
        payload.filters.typeQS ? "QS" : null,
      ]
        .filter((val) => val !== null)
        .join(",");
      if (types === "") {
        commit("emptyFilteredExperiments");
        return;
      }
      queryElements.push("experiment_type:" + types);

      // fracture (true|false)
      // fracture_mode (All modes|Mode I|Mode II|Mode III|Combined)
      if (payload.filters.withFracture && !payload.filters.withoutFracture) {
        queryElements.push("fracture:1");
        if (payload.filters.fractureMode !== null) {
          queryElements.push(`fracture_mode:${payload.filters.fractureMode}`);
        }
      } else if (
        !payload.filters.withFracture &&
        payload.filters.withoutFracture
      ) {
        queryElements.push("fracture:0");
      } else if (
        !payload.filters.withFracture &&
        !payload.filters.withoutFracture
      ) {
        commit("emptyFilteredExperiments");
        return;
      }

      // material_type_fiber_material
      if (payload.filters.fiberMaterial !== null) {
        queryElements.push(
          `material_type_fiber_material:${payload.filters.fiberMaterial}`
        );
      }

      // material_type_resin
      if (payload.filters.resin !== null) {
        queryElements.push(`material_type_resin:${payload.filters.resin}`);
      }

      // laminates_and_assemblies_stacking_sequence
      if (payload.filters.stackingSequence !== null) {
        queryElements.push(
          `laminates_and_assemblies_stacking_sequence:${payload.filters.stackingSequence}`
        );
      }

      const opts = {
        page: state.filteredExperiments.pagination.page,
        size: state.filteredExperiments.pagination.size,
        query: queryElements.join(";"),
        textSearch: payload.filters.textSearch,
      };
      this._vm.$experimentsApi.getExperimentsExperimentsGet(
        opts,
        (error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeFilteredExperiments", data);
          }
        }
      );
    },
    fetchAllFiltersValues({ commit }) {
      this._vm.$experimentsApi.getFieldDistinctExperimentsFieldDistinctGet(
        "fracture_mode",
        (error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeAllFractureMode", data);
          }
        }
      );
      this._vm.$experimentsApi.getFieldDistinctExperimentsFieldDistinctGet(
        "material_type_fiber_material",
        (error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeAllMaterialTypeFiberMaterial", data);
          }
        }
      );
      this._vm.$experimentsApi.getFieldDistinctExperimentsFieldDistinctGet(
        "material_type_resin",
        (error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeAllMaterialTypeResin", data);
          }
        }
      );
      this._vm.$experimentsApi.getFieldDistinctExperimentsFieldDistinctGet(
        "laminates_and_assemblies_stacking_sequence",
        (error, data) => {
          if (error) {
            console.error(error);
          } else {
            commit("storeAllLaminatesAndAssembliesStackingSequence", data);
          }
        }
      );
    },
  },
};
