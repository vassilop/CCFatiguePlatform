<template>
  <v-container>
    <v-card elevation="0">
      <experiment-specifications :experiment="experiment.experiment" />

      <v-data-table
        v-model="testsSelected"
        :headers="headers"
        :items="experiment.tests"
        :options.sync="options"
        :server-items-length="experiment.pagination.total"
        :loading="experiment.loadingTests"
        :footer-props="{ 'items-per-page-options': [5, 10, 15, 20, 40] }"
        item-key="id"
        @click:row="rowClick"
      >
        <template v-slot:no-data>No test for this experiment</template>
      </v-data-table>

      <v-container>
        <v-row justify="end">
          <v-btn class="ma-2" @click="goBack"> Back </v-btn>
          <v-btn
            class="ma-2"
            @click="goToTestsDashboard"
            :disabled="this.testsSelected.length === 0"
          >
            View tests Dashboard
          </v-btn>
        </v-row>
      </v-container>
    </v-card>
  </v-container>
</template>

<script>
import { mapState } from "vuex";
import ExperimentSpecifications from "@/components/ExperimentSpecifications.vue";
export default {
  name: "TestsSelection",
  components: {
    ExperimentSpecifications,
  },
  props: {
    experimentId: Number,
  },
  data() {
    return {
      testsSelected: [],
      options: {
        page: 1,
        itemsPerPage: 10,
      },

      headers: [
        { text: "Specimen Number", value: "specimen_number" },
        { text: "Stress Ratio", value: "stress_ratio" },
        { text: "Maximum Stress", value: "maximum_stress" },
        { text: "Loading Rate", value: "loading_rate" },
        { text: "Run Out", value: "run_out" },
      ],
    };
  },
  computed: {
    ...mapState("experiments", {
      experiment: "oneExperiment",
    }),
  },
  watch: {
    options: {
      handler() {
        this.fetchOneExperimentWithTests();
      },
      deep: true,
    },
  },
  methods: {
    rowClick(_item, row) {
      row.select(!row.isSelected);
    },
    goBack() {
      this.$router.go(-1);
    },
    goToTestsDashboard() {
      this.$router.push({
        name: "TestsDashboard",
        query: {
          exp: this.experimentId,
          tests: this.testsSelected.map((item) => item.id),
        },
      });
    },
    fetchOneExperimentWithTests() {
      this.$store.dispatch("experiments/fetchOneExperimentWithTests", {
        experimentId: this.experimentId,
        pagination: {
          page: this.options.page,
          size: this.options.itemsPerPage,
        },
      });
    },
  },
  created() {
    this.fetchOneExperimentWithTests();
  },
};
</script>

<style scoped lang="scss"></style>
