<template>
    <v-container fluid>
        <v-row dense>
            <v-col v-for="fshift in fshifts" :key="fshift.name" :cols="3">
                <v-card>
                    <v-card-title
                        ><router-link
                            :to="{
                                name: 'OrganismDetails',
                                params: { id: fshift.org_id },
                            }"
                            target="_blank"
                            >{{ fshift.name }}</router-link
                        ></v-card-title
                    >
                    <fshiftSignalStructure
                        v-if="fshift.signal"
                        :signal="fshift.signal"
                    ></fshiftSignalStructure>
                </v-card>
            </v-col>
        </v-row>
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
    },
};
</script>
