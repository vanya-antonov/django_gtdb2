<template>
    <v-container fluid>
        <p id="about">
            Predicted -1 frameshifting signals in chlD genes from various
            prokaryotes. The putative poly-A slippery sites and in-frame stop
            codons are marked by the blue and red colors, respectively.
        </p>
        <v-data-iterator :items="fshifts" :items-per-page="20">
            <template v-slot:default="props">
                <v-row dense>
                    <v-col
                        v-for="fshift in props.items"
                        :key="fshift.name"
                        :cols="3"
                    >
                        <v-card>
                            <v-card-title
                                ><router-link
                                    :to="{
                                        name: 'OrganismDetails',
                                        params: { id: fshift.org_id },
                                    }"
                                    target="_blank"
                                    >{{ fshift.org_name }}</router-link
                                ></v-card-title
                            >
                            <v-lazy v-model="isActive" min-height="200">
                            <fshiftSignalStructure
                                v-if="fshift.signal"
                                :signal="fshift.signal"
                            ></fshiftSignalStructure>
                            </v-lazy>
                        </v-card>
                    </v-col>
                </v-row>
            </template>
        </v-data-iterator>
    </v-container>
</template>
<script>
import HTTP from "@/http-common";
import FshiftSignalStructure from "@/components/FshiftSignalStructure";

export default {
    components: {
        FshiftSignalStructure,
    },
    data: function() {
        return {
            fshifts: [],
        };
    },
    created: function() {
        HTTP.get(`frameshifts-signals`)
            .then((response) => {
                // JSON responses are automatically parsed.
                console.log(response.data);
                this.fshifts = response.data.sort(function(a, b) {
                    return a.signal.energy - b.signal.energy;
                });
            })
            .catch((error) => {
                console.warn(error);
            });
        this.$nextTick(() => {
            window.scrollTo(0, 0);
            window.document.getElementById("about").scrollIntoView();
        });
    },
};
</script>
