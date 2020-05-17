import axios from "axios";

const HTTP = axios.create({
	baseURL: `http://localhost:8001/chelatase/api`,
});

export default HTTP;
