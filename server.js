// server.js
const express = require('express');
const bodyParser = require('body-parser');
const cors = require('cors');
const { zTest1Sample, tTest1Sample, chi2Gof, fTestTwoSample } = require('./stats');
const app = express();

app.use(cors());
app.use(bodyParser.json());
app.use(express.static('public'));

// 1) z-test endpoint
app.post('/api/ztest', (req, res) => {
  const { data, mu0, sigma, alternative } = req.body;
  const nums = data.map(Number);
  res.json(zTest1Sample(nums, Number(mu0), Number(sigma), alternative));
});

// 2) t-test endpoint
app.post('/api/ttest', (req, res) => {
  const { data, mu0, alternative } = req.body;
  const nums = data.map(Number);
  res.json(tTest1Sample(nums, Number(mu0), alternative));
});

// 3) χ²-GOF endpoint
app.post('/api/chi2', (req, res) => {
  const { observed, expected } = req.body;
  const obs = observed.map(Number), exp = expected.map(Number);
  res.json(chi2Gof(obs, exp));
});

// 4) F-test endpoint
app.post('/api/ftest', (req, res) => {
//   const data1 = req.body.data1.map(Number);
//   const data2 = req.body.data2.map(Number);
//   const alternative = req.body.alternative;
//   res.json(fTestTwoSample(data1, data2, alternative));
  // **1) Log incoming payload**
  console.log('[/api/ftest] req.body =', JSON.stringify(req.body));

  // 2) Parse arrays
  const data1 = Array.isArray(req.body.data1)
    ? req.body.data1.map(Number)
    : [];
  const data2 = Array.isArray(req.body.data2)
    ? req.body.data2.map(Number)
    : [];
  const alternative = req.body.alternative;

  // 3) Compute F-test
  const result = fTestTwoSample(data1, data2, alternative);

  // **4) Log the computed result**
  console.log('[/api/ftest] result =', result);

  // 5) Send it back
  res.json(result);
});

const PORT = process.env.PORT || 3000;
app.listen(PORT, () => console.log(`Listening on http://localhost:${PORT}`));
