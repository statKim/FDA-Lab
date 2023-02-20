ssh -p 2230 -l statKim data.opensky-network.org

show tables;

describe state_vectors_data4;

describe flights_data4;


SELECT
count(*)
FROM flights_data4
WHERE day >= 1546300800 and day <= 1548892800 and estdepartureairport = 'RKSI';