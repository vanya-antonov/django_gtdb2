<template>
    <div ref="fornaContainer" class="d-inline-flex pa-2"></div>
</template>
<script>
import { FornaContainer as PreFornaContainer } from "fornac";

export default {
    name: "fornac",
    props: {
        structure: String,
        sequence: String,
        colors: String,
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
        }

    },
    data: function() {
        return {
            fornaContainer: [],
            options: {
                structure: this.structure,
                sequence: this.sequence,
                labelInterval: 10,
                transitionDuration: 500,
                friction: 30,
            },
            colorScheme: "custom",
            customColorsText: this.colors,
        };
    },

    mounted: function() {
        this.fornaContainer[0] = new PreFornaContainer(
            this.$refs.fornaContainer,
            {
                initialSize: this.size,
                allowPanningAndZooming: this.allowPanningAndZooming,
                applyForce: this.applyForce,
            }
        );
        this.fornaContainer[0].addRNA(this.options.structure, this.options);
        this.fornaContainer[0].changeColorScheme(this.colorScheme);
        this.fornaContainer[0].addCustomColorsText(this.customColorsText);
        // this.fornaContainer[0].changeColorScheme(this.colorScheme);
        // this.fornaContainer[0].center_view();
    },
};
</script>
