<style>
.y2tick {
  display: none;
  width: 0;
}
</style>

<template>
  <div>
  <v-btn @click="plotF">plotF</v-btn>
  <div id="hideyticks">
    <vue-plotly
      ref="plotly"
      :v-if="organisms.length > 0"
      :data="heatmapData"
      :layout="layout"
      :options="options"
      @relayout="relayout"
    />
  </div>
</div>
</template>
<script>
import HTTP from "@/http-common";
import VuePlotly from "@statnett/vue-plotly";
import sendRenderEvent from "@/render-event";

function groupCategories(categories) {
  let current_category = categories[0];
  let num_current_category = 1;
  let categories_groups = [];
  let cat;
  console.log(categories)
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
  components: {
    VuePlotly,
  },

  methods: {
    plotF() {
      console.log("Run mounted")
      var update = {
          'yaxis.range': [0, this.organisms.length],   // updates the xaxis range
      };
      this.$refs.plotly.relayout(update)
    },
    afterplot(){
      console.log("afterplot")
    },
    relayout(eventData){
      console.log('relayout')
      console.log(eventData)
      if ("yaxis2.autorange" in eventData && eventData["yaxis2.autorange"]){
        var update = {
          'yaxis.range': [0, this.organisms.length],   // updates the xaxis range
        };
        this.$refs.plotly.relayout(update)
      }
    },

  },
  data: function() {
    return {
      options: {},
      organisms: [],
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
      let categories_groups = groupCategories(y_0);
      console.log(categories_groups)
      let data = categories_groups.map((cat_num) => {
        return {
          type: "bar",
          name: cat_num[0],
          x: [[''],["Taxonomy"]],
          y: [cat_num[1]],
          yaxis: "y",
        };
      });

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
        yaxis: 'y2'
        // xaxis: 'x1'
      });
      data.push({
        type: "heatmap",
        colorscale: "YlOrRd",
        reversescale: true,
        x: this.evaluesX,
        y: y,
        z: evalues_z,
        yaxis: 'y2'

        // xaxis: 'x'
      });
      console.log(data)
      return data;
    },

    layout() {
      return {
        barmode: "stack",
        showlegend:false,

        xaxis: {
          type: "multicategory",
          categoryorder: "array",
          automargin: true,
        },
        yaxis: {
          type: "linear",
          range: [0, this.organisms.length],

        },
        yaxis2: {
          type: "multicategory",
          showdividers: false,
          tickson: "boundaries",
          overlaying: 'y',
          categoryorder: "array",
          side: "left",
        },
        width: 1300,
        height: 1400,
      };
    },
  },

  created() {
    HTTP.get(`organisms`)
      .then((response) => {
        // JSON responses are automatically parsed.
        this.organisms = response.data;
        console.log("async get organisms");
        this.$nextTick(() => {
          sendRenderEvent();
        });
      })
      .catch((error) => {
        console.warn(error);
      });
  },
};
</script>
