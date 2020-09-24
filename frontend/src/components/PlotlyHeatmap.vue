<style>
.js-plotly-plot .plotly, .js-plotly-plot .plotly div {
  margin: auto !important;
}
</style>
<template>
    <vue-plotly
      ref="plotly"
      :v-if="organisms.length > 0"
      :data="heatmapData"
      :layout="layout"
      :options="options"
      @relayout="relayout"
      @click="singleClickHandler"
      @doubleclick="doubleClickHandler"
      :autoResize="false"
      class="mx-auto"
    />
</template>
<script>
import VuePlotly from "@statnett/vue-plotly";
import * as d3 from "d3";

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
      default: 700,
    },
    width: {
      type: Number,
      default: 1000,
    },
    drawBar: {
      type: Boolean,
      default: true,
    },
    removeY2ticks: {
      type: Boolean,
      default: true,
    },
    addLinksToOrganism: {
      type: Boolean,
      default: false,
    }
  },
  components: {
    VuePlotly,
  },

  methods: {
    toremoveY2ticks() {
      if (!this.removeY2ticks) {
        return;
      }
      let plot = (plot = d3.select(this.$refs.plotly.$refs.container));
      plot.attr("x", 288);
      let y2ticks = plot.selectAll(".y2tick");
      let x = y2ticks.select("text").attr("x");
      y2ticks.remove();
      plot
        .selectAll("g.y2tick2")
        .select("text")
        .attr("x", x);
    },
    relayout(eventData) {
      console.log("relayout");

      if ("yaxis2.autorange" in eventData && eventData["yaxis2.autorange"]) {
        var update = {
          "yaxis.range": [0, this.organisms.length], // updates the xaxis range
        };
        this.$refs.plotly.relayout(update);
      }
      this.toremoveY2ticks();
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
      options: {
        modeBarButtonsToRemove: [
          "zoomIn2d",
          "zoomOut2d",
          "resetScale2d",
          "pan2d",
          "select2d",
          "lasso2d",
          "hoverCompareCartesian",
          "hoverClosestCartesian",
        ],
      },
      evaluesX: [
        [
          // '',
          "Magnesium <br> chelatase <br> genes",
          "Magnesium <br> chelatase <br> genes",
          "Magnesium <br> chelatase <br> genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Chlorophyll biosynthesis genes",
          "Cobalt <br> chelatase <br> genes",
          "Cobalt <br> chelatase <br> genes",
          "Cobalt <br> chelatase <br> genes",
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
      let y_1;
      if (this.removeY2ticks) {
        y_1 = this.organisms.map((org) => org.id)
      } else if (!this.addLinksToOrganism){
        y_1 = this.organisms.map((org) => org.name);
      } else {
        y_1 = this.organisms.map((org) => {
          return `<a href="organism/`+org.id+`">`+org.name+`</a>`
        })
      }
      let y = [y_0, y_1];
      let evalues_z = this.organisms.map((org) => {
        return this.evaluesX[1].map((gene) => org.evalues[gene]);
      });
      let evalues_text = this.organisms.map((org) =>
        this.evaluesX[1].map((gene) => {
          return `
${org.name} (${org.heatmap_taxa})<br>
query protein: ${gene}<br>
E-value: 1E-${org.evalues[gene]}`;
        })
      );
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
      let ygap = y[0].length < 10? 1:0;
      console.log(ygap)
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
        ygap: ygap,
        // xaxis: 'x1'
      });
      data.push({
        type: "heatmap",
        // colorscale: "YlOrRd",
        colorscale: [
        ['0.0','red'],
        ['0.94','yellow'],
        ["1.0",'white'],
        ],
        reversescale: true,
        x: this.evaluesX,
        y: y,
        z: evalues_z,
        yaxis: "y2",
        zmax: 100,
        zmin: 0,
        ygap: ygap,
        hovertext: evalues_text,
        hoverinfo: "text",
        colorbar: {
          tickprefix: "1e-",
          showtickprefix: "all",
        },

        // xaxis: 'x'
      });
      console.log(data);
      console.log(this.layout);
      return data;
    },

    layout() {
      let automargin= !this.removeY2ticks;
      return {
        barmode: "stack",
        showlegend: false,

        xaxis: {
          type: "multicategory",
          categoryorder: "array",
          // automargin: true,
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
          automargin: automargin,
          side: "left",
          tickfont: {
            size: 15,
          }
        },
        width: this.width,
        height: this.height,
        margin: {
    l: 110,
    r: 110,
    b: 150,
    t: 10,
    // pad: 4
  }
      };
    },
  },
  mounted() {
    this.toremoveY2ticks();
  },
};
</script>
