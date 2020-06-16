<template>
    <fornac
        :structure="structure"
        :sequence="sequence"
        :colors="customColour"
        :allowPanningAndZooming="allowPanningAndZooming"
        :applyForce="applyForce"
    ></fornac>
</template>
<script>
import Fornac from "@/components/Fornac";
export default {
    name: "FshiftSignalStructre",
    props: {signal: Object,
        size: {
            type: Array,
            default: function() {
                return [300, 300];
            },
        },
        allowPanningAndZooming: {
            type: Boolean,
            default: false,
        },
        applyForce: {
            type: Boolean,
            default: false,
        }},
    components: { Fornac },
    data: function() {
        return {
            structure: this.signal.struct,
            sequence: this.signal.seq,
            stopCodonColour: "OrangeRed",
            polyAColour: "DodgerBlue"
        };
    },
    computed: {
        customColour() {
            const stop_codon_coord = this.signal.stop_codon_coord;
            let customColour = ``;
            if (stop_codon_coord) {
                customColour += `${stop_codon_coord[0] + 1}-${
                    stop_codon_coord[1]
                }:${this.stopCodonColour}`;
            }
            const poly_a_coord = this.signal.poly_a_coord;
            if (poly_a_coord) {
                customColour += ` ${poly_a_coord[0] + 1}-${
                    poly_a_coord[1]
                }:${this.polyAColour}`;
            }
            return customColour;
        },
    },
};
</script>
