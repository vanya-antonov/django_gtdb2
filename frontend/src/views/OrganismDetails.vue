<style>
.v-data-table-header th {
    white-space: nowrap;
}
</style>

<template>
    <v-col class="pt-0">
        <v-breadcrumbs
            v-if="organism.prm_dict.taxonomy"
            :items="taxonomyBreadcrumbs"
            divider="/"
            class="pa-0 ma-0"
        ></v-breadcrumbs>

        <h1>{{ organism.name }}</h1>

        <v-expansion-panels
            :value="Array.from(Array(4).keys())"
            multiple
            accordion
        >
            <v-expansion-panel v-if="organism.prm_dict.kegg_org_code">
                <v-expansion-panel-header
                    ><h3>External Links</h3></v-expansion-panel-header
                >
                <v-expansion-panel-content>
                    <ul>
                        <li>
                            <a
                                :href="
                                    `https://www.genome.jp/kegg-bin/show_pathway?${organism.prm_dict.kegg_org_code}00860`
                                "
                                >KEGG Porphyrin and chlorophyll metabolism
                                pathway</a
                            >
                        </li>
                    </ul>
                </v-expansion-panel-content>
            </v-expansion-panel>

            <v-expansion-panel>
                <v-expansion-panel-header
                    ><h3>
                        Magnesium / Cobalt chelatase subunits
                    </h3></v-expansion-panel-header
                >
                <v-expansion-panel-content>
                    <v-data-table
                        :headers="subunitsTableHeaders"
                        :items="subunitsTableItems"
                        item-key="name"
                        disable-pagination
                        fix-header
                        hide-default-footer
                        calculate-widths
                    >
                        <template v-slot:item.location="{ item }">
                            <a :href="item.hrefNcbiGeneTrack" target="_blank">{{
                                item.location
                            }}</a>
                        </template>
                    </v-data-table>
                </v-expansion-panel-content>
            </v-expansion-panel>

            <v-expansion-panel>
                <v-expansion-panel-header>
                    <h3 class="pt-3">Biosynthesis pathway genes</h3>
                </v-expansion-panel-header>
                <v-expansion-panel-content>
                    <v-row class="d-flex">
                        <v-card>
                            <v-data-table
                                :headers="pathwayGenesCountHeaders"
                                :items="b12PathwayGenesCountItems"
                                item-key="name"
                                disable-pagination
                                fix-header
                                hide-default-footer
                            >
                                <template v-slot:item.count="{ item }">
                                    <v-chip
                                        :color="getColor(item.count, 'b12')"
                                    >{{item.count}}</v-chip>
                                </template>
                            </v-data-table>
                        </v-card>
                        <v-card>
                            <v-data-table
                                :headers="pathwayGenesCountHeaders"
                                :items="chlPathwayGenesCountItems"
                                item-key="name"
                                disable-pagination
                                fix-header
                                hide-default-footer
                            >
                                <template v-slot:item.count="{ item }">
                                    <v-chip
                                        :color="getColor(item.count, 'chl')"
                                    >{{item.count}}</v-chip>
                                </template>
                            </v-data-table>
                        </v-card>
                    </v-row>
                </v-expansion-panel-content>
            </v-expansion-panel>
            <v-expansion-panel>
                <v-expansion-panel-header>
                    <h3 class="pt-3">Structure</h3>
                </v-expansion-panel-header>
                <v-expansion-panel-content>
                    <v-row><fornac></fornac></v-row>
                </v-expansion-panel-content>
            </v-expansion-panel>
        </v-expansion-panels>
    </v-col>
</template>
<script>
import HTTP from "@/http-common";
import Fornac from "@/components/Fornac";

export default {
    data: function() {
        return {
            organism: {
                prm_dict: {
                    taxonomy: [],
                },
            },
            subunitsTableHeaders: [
                { text: "Location", value: "location", width: "1%" },
                { text: "S", value: "s", width: "1%" },
                { text: "Pathway", value: "pathway", width: "1%" },
                { text: "Name", value: "name", width: "1%" },
                { text: "Len", value: "len", width: "1%" },
                { text: "FS", value: "fs", width: "1%", sortable: false },
                { text: "Gene", value: "gene", width: "1%" },
                { text: "Group", value: "group", width: "1%" },
                {
                    text: "E-value",
                    value: "evalue",
                    width: "1%",
                    sortable: true,
                },
                {
                    text: "Description",
                    value: "description",
                    width: "15%",
                    sortable: false,
                },
            ],
            pathwayGenesCountHeaders: [
                { text: "Name", value: "name", width: "1%" },
                { text: "#Genes", value: "count", width: "1%" },
            ],
        };
    },
    components: {
        Fornac,
    },
    computed: {
        taxonomyBreadcrumbs() {
            return this.organism.prm_dict.taxonomy.map((tax) => {
                return { text: tax };
            });
        },
        subunitsTableItems() {
            return this.organism.chel_subunit_feat_set.map((feat) => {
                return {
                    location: feat.prm_dict.location_str,
                    s: feat.prm_dict.chel_subunit,
                    pathway: feat.prm_dict.chel_pathway,
                    name: feat.name,
                    len: feat.prm_dict.translation_len,
                    fs: feat.frameshift_len,
                    gene: feat.prm_dict.chel_gene,
                    group: feat.prm_dict.chel_gene_group,
                    evalue: feat.prm_dict.chel_evalue,
                    description: feat.descr,
                    hrefNcbiGeneTrack: `https://www.ncbi.nlm.nih.gov/nuccore/${
                        feat.seq_id
                    }?report=graph&from=${feat.start + 1}&to=${feat.end}`,
                };
            });
        },
        b12PathwayGenesCountItems() {
            return this.organism.b12_genes_count.map((gene_count) => {
                return {
                    name: gene_count[0],
                    count: gene_count[1],
                };
            });
        },
        chlPathwayGenesCountItems() {
            return this.organism.chl_genes_count.map((gene_count) => {
                return {
                    name: gene_count[0],
                    count: gene_count[1],
                };
            });
        },
    },
    methods: {
      getColor (count, pathway) {
        var colors;
        if (pathway=="b12") {
            colors=["blue lighten-4","blue lighten-2"]
        } else {
            colors=["green lighten-4","green lighten-2"]
        }
        if (count == 1) return colors[0]
        else if (count > 1) return colors[1]
        else return "white"
      }},
    created: function() {
        HTTP.get(`organisms/` + this.$route.params.id)
            .then((response) => {
                // JSON responses are automatically parsed.
                console.log(response.data);
                this.organism = response.data;
            })
            .catch((error) => {
                console.warn(error);
            });
    },
};
</script>
