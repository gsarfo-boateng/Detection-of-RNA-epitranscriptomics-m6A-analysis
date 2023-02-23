SET GLOBAL local_infile=1;

CREATE DATABASE S2_1ug;
CREATE DATABASE S2_5ug;
-- create databases for all peak files

SHOW TABLES;

CREATE TABLE brain_s2_1ug (
  chrom VARCHAR(255),
  ChromStart INT,
  ChromEnd INT,
  strand VARCHAR(1),
  Peaklength INT,
  summit INT,
  pileup INT,
  pvlaue INT, 
  qvalue INT,
  score INT,
  PeakID VARCHAR(255)
);



USE S2_1ug ;
LOAD DATA LOCAL INFILE '/Users/georgeboateng-sarfo/Desktop/Lai_lab/m6A/miCLIP/peak_calls/TRESS_peak_call/S2/brain_s2_1ug.bed'
INTO TABLE brain_s2_1ug
FIELDS TERMINATED BY '\t' LINES TERMINATED BY '\n';


-- repeat for all peak files
