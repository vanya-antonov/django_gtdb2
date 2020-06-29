import Vue from "vue";
import VueRouter from "vue-router";
import OrganismsList from "../views/OrganismsList.vue";
import OrganismDetails from "../views/OrganismDetails.vue";
import Home from "../views/Home.vue";

// const OrganismDetails = () => import("../views/OrganismDetails.vue")

Vue.use(VueRouter);

const routes = [
  { path: "/", name: "Home", component: Home },
  {
    path: "/organisms",
    name: "OrganismsList",
    component: OrganismsList,
  },
  {
    path: "/organism/:id",
    name: "OrganismDetails",
    component: OrganismDetails,
  },
  {
    path: "/about",
    name: "About",
    // route level code-splitting
    // this generates a separate chunk (about.[hash].js) for this route
    // which is lazy-loaded when the route is visited.
    component: () =>
      import(/* webpackChunkName: "about" */ "../views/About.vue"),
  },
  {
    path: "/signals",
    name: "Signals",
    component: () =>
      import(/* webpackChunkName: "about" */ "../views/FrameShiftsSignals.vue"),
  },
];

const router = new VueRouter({
  mode: "history",
  routes,
  scrollBehavior(to, from, savedPosition) {
    if (savedPosition) {
      return savedPosition;
    } else {
      return { x: 0, y: 0 };
    }
  },
});

export default router;
