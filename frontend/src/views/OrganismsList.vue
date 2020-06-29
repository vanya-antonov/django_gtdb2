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
          ></v-autocomplete></v-col
        ><v-col md="auto">
          <v-checkbox v-model="genesToFilter.all" label="all" />
        </v-col> </v-row
      ><v-row>
        <v-autocomplete
          v-model="kingdomsToFilter"
          :items="allKingdoms"
          placeholder="Filter Kingdoms"
          chips
          multiple
        ></v-autocomplete> </v-row
      ><v-row>
        <v-autocomplete
          v-model="phylumsToFilter"
          :items="allPhylums"
          placeholder="Filter Phylums"
          chips
          multiple
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
        <v-btn @click.stop="toggleFilters">{{toggleButtonText}}</v-btn>
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
         <!-- <template v-slot:item.genotype="{ items }">
      <v-chip v-for="(item, index) in items" :key="index" color="green" dark>{{ index }}</v-chip>
    </template>
 -->      </v-data-table>
    </v-card>
  </div>
</template>
<script>
import HTTP from "@/http-common";
import sendRenderEvent from "@/render-event";

function filterGenotype(organisms, genotypes, all = false) {
  let f;
  if (all) {
    f = (org) => genotypes.every((gene) => org.genotype.includes(gene));
  } else {
    f = (org) => genotypes.some((gene) => org.genotype.includes(gene));
  }
  return organisms.filter(f);
}

function filterPhylum(organisms, phylums) {
  return organisms.filter((org) =>
    phylums.some((phylum) => org.phylum == phylum)
  );
}

function filterKingdoms(organisms, kingdoms) {
  return organisms.filter((org) =>
    kingdoms.some((kingdom) => org.kingdom == kingdom)
  );
}

export default {
  data: function() {
    return {
      drawer: null,
      search: "",
      headers: [
        { text: "Organism", value: "name" },
        { text: "Phylum", value: "phylum" },
        { text: "Kingdom", value: "kingdom" },
        { text: "Genotype", value: "genotype" },
      ],
      organisms: [],
      filteredOrganisms: [],
      genesToFilter: { values: [], all: false },
      phylumsToFilter: [],
      kingdomsToFilter: [],
      toggleButtonText: "Hide Filters",
      height: 100,
    };
  },
  computed: {
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
  },

  methods: {
    resetFilters() {
      this.genesToFilter.values = [];
      this.filteredOrganisms = this.organisms;
    },
    applyFilters() {
      let organisms = this.organisms;
      if (this.kingdomsToFilter.lenght) {
        organisms = filterKingdoms(organisms, this.kingdomsToFilter);
      }
      if (this.phylumsToFilter.lenght) {
        organisms = filterPhylum(organisms, this.phylumsToFilter);
      }
      if (this.genesToFilter.values.length) {
        organisms = filterGenotype(
          organisms,
          this.genesToFilter.values,
          this.genesToFilter.all
        );
      }

      this.filteredOrganisms = organisms;
    },
    toggleFilters() {
      this.drawer = !this.drawer
      this.toggleButtonText =  this.drawer?"Hide Filters":"Show Filters"
    }
  },

  created: function() {
    HTTP.get(`organisms`)
      .then((response) => {
        // JSON responses are automatically parsed.
        this.organisms = response.data;
        this.filteredOrganisms = this.organisms;
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
