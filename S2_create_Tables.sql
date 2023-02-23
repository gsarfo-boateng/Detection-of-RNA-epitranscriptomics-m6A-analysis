CREATE TABLE brain_s2_5ug (
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

USE S2_5ug;
LOAD DATA LOCAL infile "/Users/georgeboateng-sarfo/Desktop/Lai_lab/m6A/miCLIP/peak_calls/TRESS_peak_call/S2/brian_S2_5ug.bed"
into table brain_s2_5ug
fields terminated by "\t" LINES TERMINATED BY '\n';