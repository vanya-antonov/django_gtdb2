<style>
#textcol {
    min-width: 20rem;
}
#heatmaptext {
    min-width: 23rem;
}
</style>

<template>
    <div class="px-6">
        <v-row justify="center" class=" px-6">
            <v-col id="textcol" class="flex-shrink-1">
                <h1>Chelatase DB</h1>
                <p>
                    Magnesium chelatase chlIDH and cobalt chelatase cobNST
                    enzymes are required for biosynthesis of
                    (bacterio)chlorophyll and cobalamin (vitamin B12),
                    respectively. Each enzyme consists of large, medium and
                    small subunits. Structural and primary sequence similarities
                    indicate common evolutionary origin of the corresponding
                    subunits. This also allows prokaryotes to use taxa-specific
                    sets of genes to encode the magnesium and/or cobalt
                    chelatases.
                </p>
                <p>
                    In our recent study we analyzed 1200+ prokaryotic genomes
                    (see the <router-link to="/organisms" target="_blank">Organisms section</router-link>) and
                    identified at least 5 taxa-specific strategies that can be
                    used to synthesize these two chelatases (<a
                        href="https://pubmed.ncbi.nlm.nih.gov/32211852/"
                        >Antonov 2020</a
                    >): Interestingly, more than 100 species may utilize
                    programmed frameshifting in the chlD gene to produce both
                    the small and the medium chelatase subunits from the same
                    gene (the predicted frameshifting signals can be found in
                    the <router-link to="/signals" target="_blank">Signals section</router-link>)
                </p>
                <p>
                    Below is a graphical representation of the 5 proposed
                    strategies. For visualization purposes, a chelatase is
                    depicted as a <b>sausage</b> (small subunit) and a
                    <b>hot dog</b> (medium subunit) on a <b>plate</b> (large
                    subunit). This representation highlights the facts that each
                    of the chelatases consists of three different subunits and
                    that the small subunit is similar to a part of the medium
                    subunit (in other words, “a sausage is a part of a hot
                    dog”).
                </p>
            </v-col>
            <v-col class="flex-shrink-1 px-0">
                <chelTable></chelTable>
            </v-col>
        </v-row>
        <v-row justify="center" class=" px-6">
            <v-col id="heatmaptext">
                <h3>Prediction of the possible phenotypes</h3><ul>

<li>To determine whether the analyzed prokaryotes were able to synthesize (bacterio)chlorophyll or cobalamin (or both), we searched these genomes for genes from these biosynthesis pathways.
</li><li>On the heatmap below each row is a prokaryotic genome (1100+ species in total) and the columns correspond to presence of various genes in this genome. The species from each taxonomic group were ordered according to the 16S rRNA phylogenetic tree.
</li><li>The color of the cells corresponds to the tBLASTn -log10(E-value) produced by a set reference proteins that were used as queries. While colors indicate the gene absence. 
</li><li>All these genomes contain at least one chlD gene (see the chlD_bchD column)
</li><li>The Frameshift column indicates the presence of the frameshift mutation in the chlD gene.
</li><li>The heatmap is interactive, so the users can use mouse to determine the species name corresponding to each row as well as the -log10(E-values) for each of the identified gene in the corresponding prokaryotic genome. By clicking on each row a zoomed in version of the heatmap will appear below.</li></ul>

            </v-col>
            <v-col class="mx-auto" justify="center">
                <plotlyHeatmap
                    v-if="organisms.length"
                    :click-callback="clickCallback"
                    :organisms="organisms"
                    :height="700"
                ></plotlyHeatmap>
            </v-col>
                <v-col v-if="clickedIndex">
                <h2 align='center'>
                    <router-link
                        v-if="clickedIndex"
                        :to="{
                            name: 'OrganismDetails',
                            params: { id: organisms[clickedIndex].id },
                        }"
                        target="_blank"
                        >{{ organisms[clickedIndex].name }}</router-link
                    >
                </h2>
                <v-spacing></v-spacing>
                <!-- <v-row id="nohideticks" class="d-flex justify-center"> -->
                <plotlyHeatmap
                    v-if="clickedIndex"
                    :organisms="zoomedOrganisms"
                    :drawBar="true"
                    :click-callback="smallClickCallback"
                    :removeY2ticks="false"
                    :addLinksToOrganism="true"
                    :width="1100"
                    :height="300"
                ></plotlyHeatmap>
            </v-col>
        </v-row>
    </div>
</template>
<script>
import HTTP from "@/http-common";
import PlotlyHeatmap from "@/components/PlotlyHeatmap";
import ChelTable from "@/components/ChelTable";

export default {
    components: { PlotlyHeatmap, ChelTable },
    data: function() {
        return { organisms: [], clickedIndex: null, zoomedOrganisms: [] };
    },
    methods: {
        clickCallback(eventData) {
            console.log("callback");
            console.log(eventData);
            if (eventData.points[0].data.type == "heatmap") {
                this.clickedIndex = eventData.points[0].pointIndex[0];
                this.setZoomedOrganisms();
            }
        },
        smallClickCallback(eventData) {
            console.log(eventData);
            // if (eventData.points[0].data.type == "heatmap") {
            //  this.clickedIndex = eventData.points[0].pointIndex[0];
            // }
        },
        setZoomedOrganisms() {
            this.zoomedOrganisms = this.organisms.slice(
                Math.max(0, this.clickedIndex - 3),
                Math.min(this.clickedIndex + 4, this.organisms.length)
            );
            console.log(this.zoomedOrganisms);
        },
    },

    created() {
        HTTP.get(`organisms`)
            .then((response) => {
                // JSON responses are automatically parsed.
                this.organisms = response.data;
            })
            .catch((error) => {
                console.warn(error);
            });
    },
};
</script>
