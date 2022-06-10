-- MySQL dump 10.14  Distrib 5.5.56-MariaDB, for Linux (x86_64)
--
-- Host: mysql-pfam-rel    Database: pfam_35_0
-- ------------------------------------------------------
-- Server version	5.6.36-log

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `pdb_residue_data`
--

DROP TABLE IF EXISTS `pdb_residue_data`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!40101 SET character_set_client = utf8 */;
CREATE TABLE `pdb_residue_data` (
  `pdb_id` varchar(5) NOT NULL,
  `chain` varchar(4) CHARACTER SET latin1 COLLATE latin1_general_cs DEFAULT NULL,
  `serial` int(10) DEFAULT NULL,
  `pdb_res` char(3) DEFAULT NULL,
  `pdb_seq_number` int(10) DEFAULT NULL,
  `pdb_insert_code` varchar(1) DEFAULT NULL,
  `observed` int(1) DEFAULT NULL,
  `dssp_code` varchar(4) DEFAULT NULL,
  `pfamseq_acc` varchar(10) NOT NULL,
  `pfamseq_res` char(3) DEFAULT NULL,
  `pfamseq_seq_number` int(10) DEFAULT NULL,
  KEY `pdb_id` (`pdb_id`),
  KEY `serial_idx` (`serial`),
  KEY `obs_idx` (`observed`),
  KEY `pfamseq_seq_number_idx` (`pfamseq_seq_number`),
  KEY `chain_idx` (`chain`),
  KEY `fk_pdb_residue_data_uniprot1_idx` (`pfamseq_acc`),
  CONSTRAINT `FK_pdb_residue_data_1` FOREIGN KEY (`pdb_id`) REFERENCES `pdb` (`pdb_id`) ON DELETE CASCADE ON UPDATE NO ACTION,
  CONSTRAINT `fk_pdb_residue_data_uniprot1` FOREIGN KEY (`pfamseq_acc`) REFERENCES `uniprot` (`uniprot_acc`) ON DELETE CASCADE ON UPDATE NO ACTION
) ENGINE=MyISAM DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2021-11-10 14:28:18
