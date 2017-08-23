/* gthDB_make.sql
 * Michael E Sparks (mespar1@iastate.edu)
 *
 * Last Modified: 22 October 2005
 * 226 lines, 7732 characters
 *
 * This is a script of SQL commands to create gthDB,
 * a database for warehousing GSQ/GTH output's core data.
 *
 * Copyright (c) 2003,2004,2005 Michael E Sparks
 * All rights reserved.

 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:

     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright
       notice, this list of conditions and the following disclaimer in
       the documentation and/or other materials provided with the
       distribution.
     * Neither the name of the Iowa State University nor the names of
       its contributors may be used to endorse or promote products
       derived from this software without specific prior written permission.

 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* Create the database */
CREATE DATABASE gthDB;

/* Now generate the tables... */

CREATE TABLE gthDB.header (
  /* Corresponds to each case of [GSQ|GTH] output */
  run SMALLINT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  /* Was it GeneSeqer, GenomeThreader, or other that produced the data? */
  program VARCHAR(25),
  /* Version control */
  version VARCHAR(25),
  /* When was the data produced? */
  daterun VARCHAR(50),
  /* What species-specific splice site model was used in the alignment? */
  splicemodel VARCHAR(50),
  /* Run-time parameters used */
  parms TEXT,
  /* Source file(s) for template sequences */
  template_files TEXT,
  /* Source file(s) for reference sequences */
  ref_files TEXT,
  /* Was the reference ESTcDNA, Protein, or Mixed? */
  ref_type ENUM("ESTcDNA","Protein","Mixed")
);

CREATE TABLE gthDB.alignments (
  /* ties with header.run */
  run SMALLINT UNSIGNED NOT NULL,
  /* unique [ESTcDNA|Protein]_information entries */
  alignID MEDIUMINT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  /* Name of reference sequence */
  refID VARCHAR(40),
  /* Orientation of reference sequence */
  refstrand ENUM("+","-"),
  /* Source file of reference sequence */
  ref_file VARCHAR(100),
  /* Name of template sequence */
  templateID VARCHAR(40),
  /* Orientation of template sequence */
  templatestrand ENUM("+","-"),
  /* Source file of template sequence */
  template_file VARCHAR(100),
  /* Start of template sequence used in the alignment */
  temp_start INT UNSIGNED NOT NULL,
  /* Stop of template sequence used in the alignment */
  temp_stop INT UNSIGNED NOT NULL,
  /* Score of the alignment */
  total_score FLOAT,
  /* Total length of exons from the alignment */
  length_exons MEDIUMINT UNSIGNED,
  /* Total count of exons in the alignment */
  tot_exons TINYINT UNSIGNED NOT NULL,
  /* Percentage of coverage */
  coverage FLOAT,
  /* Which item was covered the most? */
  type_high_cov ENUM("G","C","P"),
  /* Was a high quality PGS reported for the alignment? */
  hiqual_bool ENUM("FALSE","TRUE")
);

CREATE TABLE gthDB.align_exons (
  /* Ties with alignments.alignID */
  alignID MEDIUMINT UNSIGNED NOT NULL,
  /* Serial number of exon in the alignment */
  serial TINYINT UNSIGNED NOT NULL,
  /* Start of exon region of template */
  gen_start INT UNSIGNED NOT NULL,
  /* Stop of exon region of template */
  gen_stop INT UNSIGNED NOT NULL,
  /* Start of exon region of reference sequence */
  ref_start SMALLINT UNSIGNED NOT NULL,
  /* Stop of exon region of reference sequence */
  ref_stop SMALLINT UNSIGNED NOT NULL,
  /* Score of the exon alignment */
  ref_score FLOAT
);

CREATE TABLE gthDB.align_introns (
  /* Ties with alignments.alignID */
  alignID MEDIUMINT UNSIGNED NOT NULL,
  /* Serial number of intron in the alignment */
  serial TINYINT UNSIGNED NOT NULL,
  /* Start of intron region of template */
  int_start INT UNSIGNED NOT NULL,
  /* Stop of intron region of template */
  int_stop INT UNSIGNED NOT NULL,
  /* Probability of donor site */
  don_prob FLOAT,
  /* Score of donor site */
  don_score FLOAT,
  /* Probability of acceptor site */
  acc_prob FLOAT,
  /* Score of acceptor site */
  acc_score FLOAT
);

CREATE TABLE gthDB.hiqual_align_exons (
  /* Ties with alignments.alignID */
  alignID MEDIUMINT UNSIGNED NOT NULL,
  /* Serial number of high quality exon in the alignment */
  serial TINYINT UNSIGNED NOT NULL,
  /* Start of exon region of template */
  gen_start INT UNSIGNED NOT NULL,
  /* Stop of exon region of template */
  gen_stop INT UNSIGNED NOT NULL
);

CREATE TABLE gthDB.predgenloc (
  /* Ties with header.run */
  run SMALLINT UNSIGNED NOT NULL,
  /* Unique predicted gene location entries */
  pglID MEDIUMINT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  /* Start of predicted gene location on template sequence */
  pgl_start INT UNSIGNED NOT NULL,
  /* Stop of predicted gene location on template sequence */
  pgl_stop INT UNSIGNED NOT NULL
);

CREATE TABLE gthDB.altgenstruct (
  /* Ties with predgenloc.pglID */
  pglID MEDIUMINT UNSIGNED NOT NULL,
  /* Unique AGS entries */
  agsID MEDIUMINT UNSIGNED NOT NULL PRIMARY KEY AUTO_INCREMENT,
  /* Supporting ESTcDNA/Protein evidence */
  support VARCHAR(250),
  /* Total count of exons in the AGS */
  tot_exons TINYINT UNSIGNED NOT NULL
);

CREATE TABLE gthDB.ags_exons (
  /* Ties with altgenstruct.agsID */
  agsID MEDIUMINT UNSIGNED NOT NULL,
  /* Serial number of exon in putative gene structure */
  serial TINYINT UNSIGNED NOT NULL,
  /* Start of putative exon on template */
  e_start INT UNSIGNED NOT NULL,
  /* Stop of putative exon on template */
  e_stop INT UNSIGNED NOT NULL,
  /* Score of the exon */
  e_score FLOAT
);

CREATE TABLE gthDB.ags_introns (
  /* Ties with altgenstruct.agsID */
  agsID MEDIUMINT UNSIGNED NOT NULL,
  /* Serial number of intron in putative gene structure */
  serial TINYINT UNSIGNED NOT NULL,
  /* Start of putative intron on template */
  i_start INT UNSIGNED NOT NULL,
  /* Stop of putative intron on template */
  i_stop INT UNSIGNED NOT NULL,
  /* Probability of putative donor site */
  don_prob FLOAT,
  /* Probability of putative acceptor site */
  acc_prob FLOAT
);

CREATE TABLE gthDB.prob_orfs (
  /* Ties with altgenstruct.agsID */
  agsID MEDIUMINT UNSIGNED NOT NULL,
  /* Name of template sequence */
  templateID VARCHAR(40),
  /* Orientation of reference sequence */
  template_strand ENUM("+","-"),
  /* Exon coordinates */
  exon_coors TEXT,
  /* Frame of the pps */ 
  frame TINYINT UNSIGNED NOT NULL,
  /* Length of predicted polypeptide open reading frame */
  pep_length MEDIUMINT UNSIGNED NOT NULL,
  /* Sequence of the predicted polypeptide */
  pep_seq TEXT
);

/* Handy statements for resetting gthDB
use gthDB;
truncate header;
truncate alignments;
truncate align_exons;
truncate align_introns;
truncate hiqual_align_exons;
truncate predgenloc;
truncate altgenstruct;
truncate ags_exons;
truncate ags_introns;
truncate prob_orfs;
*/
