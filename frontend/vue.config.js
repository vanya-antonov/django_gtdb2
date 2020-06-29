module.exports = {
  transpileDependencies: ["vuetify"],
  publicPath: process.env.NODE_ENV === 'development' ? '/' : '/chelatase/',
  configureWebpack: {
    module: {
      rules: [
        {
          test: /\.js$/,
          use: [
            "ify-loader",
            "transform-loader?plotly.js/tasks/compress_attributes.js",
          ],
        },
      ],
    },
  },
};
