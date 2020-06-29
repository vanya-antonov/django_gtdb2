<template>
  <div>
    <v-navigation-drawer
      v-model="drawer"
      app
      clipped
      right
      width="30%"
      class="pa-2"
      ><v-col>
        <v-row>
          <v-col> <v-btn @click="applyFilters">Apply Filters</v-btn></v-col>
          <v-col><v-btn @click="resetFilters">Reset Filters</v-btn></v-col>
        </v-row>
        <v-row
          ><v-col
            ><v-autocomplete
              v-model="genesToFilter.values"
              :items="allGenotypeGenes"
              placeholder="Filter Genotypes"
              chips
              multiple
              :error="!genesInItems"
            ></v-autocomplete></v-col
          ><v-col md="auto">
            <v-radio-group v-model="genesToFilter.type">
      <v-radio
        :key="0"
        label="any"
        value="any"
      ></v-radio>
      <v-radio
        :key="1"
        label="all"
        value="all"
      ></v-radio>
      <v-radio
        :key="2"
        label="exact"
        value="exact"
      ></v-radio>
    </v-radio-group>
          </v-col> </v-row
        ><v-row>
          <v-autocomplete
            v-model="kingdomsToFilter"
            :items="allKingdoms"
            placeholder="Filter Kingdoms"
            chips
            multiple
            :error="!kingdomsInItems"
          ></v-autocomplete> </v-row
        ><v-row>
          <v-autocomplete
            v-model="phylumsToFilter"
            :items="allPhylums"
            placeholder="Filter Phylums"
            chips
            multiple
            :error="!phylumsInItems"
          ></v-autocomplete>
        </v-row>
      </v-col>
    </v-navigation-drawer>

    <v-card class="mt-auto">
      <v-card-title class="mt-0 pt-0">
        <v-text-field
          width="200px"
          v-model="search"
          append-icon="mdi-magnify"
          label="Search"
          single-line
          hide-details
        ></v-text-field>
        <v-spacer></v-spacer>
        <v-btn @click.stop="toggleFilters">{{ toggleButtonText }}</v-btn>
      </v-card-title>
      <v-data-table
        :headers="headers"
        :items="filteredOrganisms"
        :search="search"
        item-key="name"
        fixed-header
        disable-pagination
        show-group-by
        hide-default-footer
        height="calc(100vh - 130px)"
        style="max-height: calc(100vh - 144px)"
      >
        <template v-slot:item.name="{ item }">
          <router-link
            :to="{
              name: 'OrganismDetails',
              params: { id: item.id },
            }"
            target="_blank"
            >{{ item.name }}</router-link
          >
        </template>
      </v-data-table>
    </v-card>
  </div>
</template>
<script>
import HTTP from "@/http-common";
import sendRenderEvent from "@/render-event";

function filterGenotype(organisms, genotypes, type) {
  console.log("filterGenotype");
  console.log(type)
  let f;
  if (type == 'all') {
    f = (org) => genotypes.every((gene) => org.genotype.includes(gene));
  } else if (type== 'any') {
    f = (org) => genotypes.some((gene) => org.genotype.includes(gene));
  } else if (type == 'exact') {
    f = (org) => (genotypes.sort().join(', ') == org.genotypeText)
  } else {
    throw "ValueError"
  }
  return organisms.filter(f);
}

function filterPhylum(organisms, phylums) {
  console.log("filterPhylum");
  return organisms.filter((org) =>
    phylums.some((phylum) => org.phylum == phylum)
  );
}

function filterKingdoms(organisms, kingdoms) {
  console.log("filterKingdoms");

  return organisms.filter((org) =>
    kingdoms.some((kingdom) => org.kingdom == kingdom)
  );
}

export default {
  props: {
    genesToFilterProp: {
      type: Object,
      default: () => ({ values: [], type: 'any' }),
    },
    phylumsToFilterProp: {
      type: Array,
      default: () => [],
    },
    kingdomsToFilterProp: {
      type: Array,
      default: () => [],
    },
  },
  data: function() {
    return {
      drawer: null,
      search: "",
      headers: [
        { text: "Organism", value: "name" },
        { text: "Phylum", value: "phylum" },
        { text: "Kingdom", value: "kingdom" },
        { text: "Genotype", value: "genotypeText" },
      ],
      organisms: [],
      filteredOrganisms: [],
      genesToFilter: this.genesToFilterProp,
      phylumsToFilter: this.phylumsToFilterProp,
      kingdomsToFilter: this.kingdomsToFilterProp,
      height: 100,
    };
  },
  computed: {
    genesInItems() {
      return this.genesToFilter.values.every((gene) =>
        this.allGenotypeGenes.includes(gene)
      );
    },
    phylumsInItems() {
      return this.phylumsToFilter.every((phylum) =>
        this.allPhylums.includes(phylum)
      );
    },
    kingdomsInItems() {
      return this.phylumsToFilter.every((phylum) =>
        this.allPhylums.includes(phylum)
      );
    },

    allGenotypeGenes() {
      return Array.from(
        new Set(this.organisms.map((org) => org.genotype).flat())
      ).sort();
    },
    allPhylums() {
      return Array.from(
        new Set(this.organisms.map((org) => org.phylum).flat())
      );
    },
    allKingdoms() {
      return Array.from(
        new Set(this.organisms.map((org) => org.kingdom).flat())
      );
    },
    toggleButtonText() {
      return this.drawer ? "Hide Filters" : "Show Filters";
    },
  },

  methods: {
    resetFilters() {
      this.genesToFilter.values = [];
      this.filteredOrganisms = this.organisms;
    },
    applyFilters() {
      let organisms = this.organisms;
      console.log([
        this.kingdomsToFilter,
        this.phylumsToFilter,
        this.genesToFilter,
      ]);
      console.log(this.kingdomsToFilter.length);
      if (this.kingdomsToFilter.length) {
        console.log("filter kingdoms");
        organisms = filterKingdoms(organisms, this.kingdomsToFilter);
      }
      if (this.phylumsToFilter.length) {
        organisms = filterPhylum(organisms, this.phylumsToFilter);
      }
      if (this.genesToFilter.values.length) {
        organisms = filterGenotype(
          organisms,
          this.genesToFilter.values,
          this.genesToFilter.type
        );
      }

      this.filteredOrganisms = organisms;
    },
    toggleFilters() {
      this.drawer = !this.drawer;
    },
  },

  created: function() {
    HTTP.get(`organisms`)
      .then((response) => {
        // JSON responses are automatically parsed.
        this.organisms = response.data.map((org) => ({
          ...org,
          genotypeText: org.genotype.sort().join(", "),
        }));
        this.filteredOrganisms = this.organisms;
        this.applyFilters();
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
