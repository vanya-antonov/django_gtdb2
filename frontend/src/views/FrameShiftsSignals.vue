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
                    <v-lazy
                        v-model="isActive"
                        min-height="200"
                        transition="fade-transition"
                    >
                        <fshiftSignalStructure
                            v-if="fshift.signal"
                            :signal="fshift.signal"
                        ></fshiftSignalStructure>
                    </v-lazy>
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
    mounted: function() {
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
        // this.fshifts = [
        //     {
        //         org_name: "Methanocaldococcus sp. FS406-22",
        //         org_id: 84,
        //         name: "NC_013887.1:1224254:-1",
        //         poly_a_slippery: [1224256, 1224265],
        //         signal: {
        //             struct:
        //                 "..........................(((((((..(........)..)))))))...................................................((((((((((((....((((((...)))))).....))))))))...))))(((((((.......))))))).",
        //             seq:
        //                 "AAAAAAAAACAUGAUGAAAUAAGAAAUGAGUUUGAAGAAGAAAAUGAGGAUUCAAAUAAUCAAAAUAAUAAUAACAACUCUAAUAACCAAAAUGAAGAUACACCCGGGGACUUUGAAAGAACGUUUGGCAUAGAUGAGAGCUUUAAAGUAAAUCCCAAGCUUAUACAAUUUAAGCUUA",
        //             poly_a_coord: [0, 9],
        //             stop_codon_coord: [19, 22],
        //             energy: -0.11741573033707864,
        //         },
        //         strand: 1,
        //     },
        //     {
        //         org_name: "Methanocaldococcus fervens AG86",
        //         org_id: 85,
        //         name: "NC_013156.1:1122957:-1",
        //         poly_a_slippery: [1122935, 1122946],
        //         signal: {
        //             struct:
        //                 "..................................................................................(((........)))........((((((((((((((...((((((...))))))...)).))))))))..)))).",
        //             seq:
        //                 "AAAAAAAAAAACAAGAAUAAAUUAAAUAAUGAAUCAAAUAAUAAUGAGAAUAAUAAUCCAAAUAAUGAACAUGAAAAUAAUAAUCAAAACCAAGAUGAAAAUACUGGAGAUUUUGAGCAAACAUUUGGUAUAGAUGAGAGCGUUAAAGUUAAUCCAA",
        //             poly_a_coord: [0, 11],
        //             stop_codon_coord: [22, 25],
        //             energy: -0.07579617834394904,
        //         },
        //         strand: 1,
        //     },
        // ];
        console.log("get fshifts");
        console.log(this.fshifts);
    },
};
</script>
