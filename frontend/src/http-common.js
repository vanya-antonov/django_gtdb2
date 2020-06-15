import axios from "axios";

const HTTP = axios.create({
	baseURL:
		`api` ||
		`http://localhost:8001/chelatase/api`,
});

export default HTTP;
