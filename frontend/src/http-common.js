import axios from "axios";

const HTTP = axios.create({
	baseURL:
		process.env.VUE_APP_API_BASE_URL ||
		`http://localhost:8001/chelatase/api`,
});

export default HTTP;
