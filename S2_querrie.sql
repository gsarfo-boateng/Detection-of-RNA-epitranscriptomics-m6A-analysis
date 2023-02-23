CREATE DATABASE unique_regions;

USE unique_regions;

CREATE TABLE regions (
  chromosome VARCHAR(50),
  start INT(10),
  end INT(10),
  strand CHAR(1),
  peak_file VARCHAR(50)
);

INSERT INTO regions
SELECT distinct chrom, chromStart, ChromEnd, strand FROM S2_1ug.brain_s2_1ug
UNION
SELECT DISTINCT chrom, chromstart, chromend, strand FROM S2_5ug.brain_s2_5ug;
-- 'S2_1ug' AS peak_file, 'S2_5ug' AS peak_file

-- repeat for all peak files

SELECT DISTINCT chrom FROM S2_5ug.brain_s2_5ug;

ALTER TABLE S2_5ug.brain_s2_5ug ADD COLUMN redundancy_regions INT DEFAULT 0;

-- update new columns with 
UPDATE S2_5ug.brain_s2_5ug
SET redundancy_count = (
  SELECT COUNT(*)
  FROM S2_5ug.brain_s2_5ug AS p2
  WHERE S2_5ug.brain_s2_5ug.chrom = p2.chrom
    AND S2_5ug.brain_s2_5ug.ChromStart = p2.ChromStart
    AND S2_5ug.brain_s2_5ug.ChromEnd = p2.ChromEnd
    AND S2_5ug.brain_s2_5ug.strand = p2.strand
    AND S2_5ug.brain_s2_5ug.PeakID = p2.peakID
);

SELECT chrom, ChromStart, ChromEnd, COUNT(*) AS redundancy_count
FROM S2_5ug.brain_s2_5ug
GROUP BY chrom, ChromStart, ChromEnd, strand, PeakID
HAVING COUNT(*) > 1;

-- Select jus the redundant regions from the 5ug 
SELECT chrom, chromStart, chromEnd, COUNT(*) AS redundancy
FROM S2_5ug.brain_s2_5ug
GROUP BY chrom, chromStart, chromEnd
HAVING COUNT(*) > 1;


-- Select just the redundant regions from the s2_1ug 
SELECT chrom, chromStart, chromEnd, COUNT(*) AS redundancy
FROM S2_1ug.brain_s2_1ug
GROUP BY chrom, chromStart, chromEnd
HAVING COUNT(*) > 1;

-- Determine the percentage of redundant genomic regions in S2_1ug
SELECT COUNT(DISTINCT concat_ws('_', chrom, chromStart, chromEnd)) AS unique_regions,
       COUNT(*) AS total_regions,
       (COUNT(*) - COUNT(DISTINCT concat_ws('_', chrom, chromStart, chromEnd))) AS redundant_regions,
       ROUND(((COUNT(*) - COUNT(DISTINCT concat_ws('_', chrom, chromStart, chromEnd))) / COUNT(*) * 100), 2) AS percentage_redundant
FROM S2_1ug.brain_s2_1ug;

-- Determine statistics for redundant regions in S2_5ug
SELECT COUNT(DISTINCT concat_ws('_', chrom, chromStart, chromEnd)) AS unique_regions,
       COUNT(*) AS total_regions,
       (COUNT(*) - COUNT(DISTINCT concat_ws('_', chrom, chromStart, chromEnd))) AS redundant_regions,
       ROUND(((COUNT(*) - COUNT(DISTINCT concat_ws('_', chrom, chromStart, chromEnd))) / COUNT(*) * 100), 2) AS percentage_redundant
FROM S2_5ug.brain_s2_5ug;

