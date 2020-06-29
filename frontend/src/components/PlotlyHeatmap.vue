<style>
.y2tick {
  display: none;
  width: 0;
}
</style>

<template>
  <div id="hideyticks">
    <vue-plotly
      ref="plotly"
      :v-if="organisms.length > 0"
      :data="heatmapData"
      :layout="layout"
      :options="options"
      @relayout="relayout"
      @click="singleClickHandler"
      @doubleclick="doubleClickHandler"
    />
  </div>
</template>
<script>
import VuePlotly from "@statnett/vue-plotly";

function groupCategories(categories) {
  let current_category = categories[0];
  let num_current_category = 1;
  let categories_groups = [];
  let cat;
  for (cat of categories.slice(1)) {
    if (cat != current_category) {
      categories_groups.push([current_category, num_current_category]);
      num_current_category = 1;
      current_category = cat;
    } else {
      num_current_category += 1;
    }
  }
  categories_groups.push([current_category, num_current_category]);
  return categories_groups;
}

export default {
  props: {
    organisms: Array,

    clickCallback: { type: Function, default: function() {} },
    height: {
      type: Number,
      default: 500,
    },
    width: {
      type: Number,
      default: 1300,
    },
    drawBar: {
      type: Boolean,
      default: true,
    },
  },
  components: {
    VuePlotly,
  },

  methods: {
    relayout(eventData) {
      console.log("relayout");

      if ("yaxis2.autorange" in eventData && eventData["yaxis2.autorange"]) {
        var update = {
          "yaxis.range": [0, this.organisms.length], // updates the xaxis range
        };
        this.$refs.plotly.relayout(update);
      }
    },
    singleClickHandler(eventData) {
      const interval = 500;
      let t0 = Date.now();

      // make sure click isn't 2nd click in a double-click
      if (t0 - this.doubleClickTime > interval) {
        // wait enough time for double-click to make sure this isn't the first click in a double-click
        setTimeout(() => {
          if (t0 - this.doubleClickTime > interval) {
            this.clickCallback(eventData);
          }
        }, interval);
      }
    },

    doubleClickHandler() {
      this.doubleClickTime = Date.now();
      console.log("double click");
    },
  },
  data: function() {
    return {
      doubleClickTime: 0,
      options: {},
      evaluesX: [
        [
          // '',
          "Magnesium chelatase genes",
          "Magnesium chelatase genes",
          "Magnesium chelatase genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Cobalt chelatase genes",
          "Cobalt chelatase genes",
          "Cobalt chelatase genes",
          "Vitamin B12 biosynthesis genes",
          "Vitamin B12 biosynthesis genes",
          "Vitamin B12 biosynthesis genes",
          "Vitamin B12 biosynthesis genes",
          "Vitamin B12 biosynthesis genes",
          "Vitamin B12 biosynthesis genes",
        ],
        [
          // "Frameshift",
          "chlH_bchH",
          "chlD_bchD",
          "chlI_bchI",
          "bchE",
          "chlB_bchB",
          "chlG_bchG",
          "chlL_bchL",
          "chlM_bchM",
          "chlN_bchN",
          "cobN",
          "cobT",
          "cobS",
          "cobD_cobC",
          "cobO",
          "cobP_cobU",
          "cobQ",
          "cobV_cobS",
          "cysG_cobA",
        ],
      ],
      fshiftsX: [[""], ["Frameshift"]],
    };
  },
  computed: {
    heatmapData() {
      let y_0 = this.organisms.map((org) => org.heatmap_taxa);
      let y_1 = this.organisms.map((org) => org.name);
      let y = [y_0, y_1];
      let evalues_z = this.organisms.map((org) => {
        return this.evaluesX[1].map((gene) => org.evalues[gene]);
      });
      let fshifts_z = this.organisms.map((org) => [org.num_fshifts]);
      let data = [];
      if (this.drawBar) {
        let categories_groups = groupCategories(y_0);

        data = categories_groups.map((cat_num) => {
          return {
            type: "bar",
            name: cat_num[0],
            x: [[""], ["Taxonomy"]],
            y: [cat_num[1]],
            yaxis: "y",
          };
        });
      }

      data.push({
        type: "heatmap",
        colorscale: [
          ["0.0", "rgb(255,255,255)"],
          ["1.0", "black"],
        ],
        x: this.fshiftsX,
        y: y,
        z: fshifts_z,
        showscale: false,
        yaxis: "y2",
        zmax: 2,
        zmin: 0,
        // xaxis: 'x1'
      });
      data.push({
        type: "heatmap",
        colorscale: "YlOrRd",
        reversescale: true,
        x: this.evaluesX,
        y: y,
        z: evalues_z,
        yaxis: "y2",
        zmax: 100,
        zmin: 0,

        // xaxis: 'x'
      });
      console.log(data);
      console.log(this.layout);
      return data;
    },

    layout() {
      return {
        barmode: "stack",
        showlegend: false,

        xaxis: {
          type: "multicategory",
          categoryorder: "array",
          automargin: true,
          mirror: true,
        },
        yaxis: {
          type: "linear",
          range: [0, this.organisms.length],
          visible: false,
        },
        yaxis2: {
          type: "multicategory",
          showdividers: false,
          tickson: "boundaries",
          overlaying: "y",
          categoryorder: "array",
          automargin: true,
          side: "left",
        },
        width: this.width,
        height: this.height,
      };
    },
  },
};
</script>
