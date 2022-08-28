// License:
// SurvivalGWAS_SV v1.3.2 (Multithread process) - SurvivalGWAS_SV is a single variant analytics tool for GWAS.
// Copyright © 2016, Hamzah Syed, University of Liverpool.
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


// Accord.NET Applications
// http://accord-framework.net
//
// Copyright © 2009-2014, César Souza
// All rights reserved. 3-BSD License:

//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are met:
//
//      * Redistributions of source code must retain the above copyright
//        notice, this list of conditions and the following disclaimer.
//
//      * Redistributions in binary form must reproduce the above copyright
//        notice, this list of conditions and the following disclaimer in the
//        documentation and/or other materials provided with the distribution.
//
//      * Neither the name of the Accord.NET Framework authors nor the
//        names of its contributors may be used to endorse or promote products
//        derived from this software without specific prior written permission.
// 
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
//  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Data;
using Accord.Math;
using Accord.Statistics.Analysis;
using Accord.Statistics.Testing;
using NDesk.Options;
using Accord.Math.Integration;
using System.IO.Compression;
using System.Threading;
using System.Threading.Tasks;

namespace SurvivalGWAS
{
    class Program
    {
        public static string cores;
        public static string file;
        public static string file2;
        public static string Stime;
        public static string cens;
        public static string int1;
        public static string covariates;
        public static string start;
        public static string stop;
        public static string print;
        public static string SaveFile;
        public static string method;
        public static string num;
        public static string num2;
        public static string chrome;
        public static bool show_help = false;
        public static int chr;
        public static string decompressedFileName;
        public static List<long> indexes;
        public static ReaderWriterLock rwl = new ReaderWriterLock();
        // public static object mylock = new object();
        [ThreadStatic]
        private static ProportionalHazardsAnalysis pha;
        [ThreadStatic]
        private static ProportionalHazardsAnalysis pha2;
        [ThreadStatic]
        static int number;
        [ThreadStatic]
        static int number2;
        [ThreadStatic]
        static DataTable dos;
        [ThreadStatic]
        static DataTable dt;
        [ThreadStatic]
        static DataTable table;
        [ThreadStatic]
        static DataTable independent;
        [ThreadStatic]
        static DataRow dr;
        [ThreadStatic]
        static String[] newnames;
        [ThreadStatic]
        static string[][] splitFileContents2;
        [ThreadStatic]
        static double[] time;
        [ThreadStatic]
        static int[] censor;
        [ThreadStatic]
        static double[][] inputer;
        [ThreadStatic]
        static List<string> names;
        [ThreadStatic]
        static List<string> names2;
        [ThreadStatic]
        static DataColumn intColumn;
        [ThreadStatic]
        static double allele;
        [ThreadStatic]
        static double infoscore;
        [ThreadStatic]
        static string[] split;
        [ThreadStatic]
        static string[] split2;
        [ThreadStatic]
        static string[] values;
        [ThreadStatic]
        static List<double> NN;
        [ThreadStatic]
        static List<double> PP;
        [ThreadStatic]
        static List<double> FF;
        [ThreadStatic]
        static List<double> FPsub;
        [ThreadStatic]
        static List<double> gsum2;
        [ThreadStatic]
        static string lrt;
        [ThreadStatic]
        static string[] a;
        [ThreadStatic]
        static int valid;
        [ThreadStatic]
        static int itercount;
        [ThreadStatic]
        static double b0;
        [ThreadStatic]
        static double b1;
        [ThreadStatic]
        static double b2;
        [ThreadStatic]
        static double b3;
        [ThreadStatic]
        static double b4;
        [ThreadStatic]
        static double b5;
        [ThreadStatic]
        static double b6;
        [ThreadStatic]
        static double b7;
        [ThreadStatic]
        static double b8;
        [ThreadStatic]
        static double b9;
        [ThreadStatic]
        static double b10;
        [ThreadStatic]
        static double b11;
        [ThreadStatic]
        static double sig;
        [ThreadStatic]
        static double sumda;
        [ThreadStatic]
        static double sumdb;
        [ThreadStatic]
        static double sumdc;
        [ThreadStatic]
        static double sumdd;
        [ThreadStatic]
        static double sumde;
        [ThreadStatic]
        static double sumdf;
        [ThreadStatic]
        static double sumdg;
        [ThreadStatic]
        static double sumdh;
        [ThreadStatic]
        static double sumdi;
        [ThreadStatic]
        static double sumdj;
        [ThreadStatic]
        static double sumdk;
        [ThreadStatic]
        static double sumdl;
        [ThreadStatic]
        static double sumdm;
        [ThreadStatic]
        static double sumd2a;
        [ThreadStatic]
        static double sumd2b;
        [ThreadStatic]
        static double sumd2c;
        [ThreadStatic]
        static double sumd2d;
        [ThreadStatic]
        static double sumd2e;
        [ThreadStatic]
        static double sumd2f;
        [ThreadStatic]
        static double sumd2g;
        [ThreadStatic]
        static double sumd2h;
        [ThreadStatic]
        static double sumd2i;
        [ThreadStatic]
        static double sumd2j;
        [ThreadStatic]
        static double sumd2k;
        [ThreadStatic]
        static double sumd2l;
        [ThreadStatic]
        static double sumd2m;
        [ThreadStatic]
        static double sumdab;
        [ThreadStatic]
        static double sumdac;
        [ThreadStatic]
        static double sumdad;
        [ThreadStatic]
        static double sumdae;
        [ThreadStatic]
        static double sumdaf;
        [ThreadStatic]
        static double sumdag;
        [ThreadStatic]
        static double sumdah;
        [ThreadStatic]
        static double sumdai;
        [ThreadStatic]
        static double sumdaj;
        [ThreadStatic]
        static double sumdak;
        [ThreadStatic]
        static double sumdal;
        [ThreadStatic]
        static double sumdbc;
        [ThreadStatic]
        static double sumdbd;
        [ThreadStatic]
        static double sumdbe;
        [ThreadStatic]
        static double sumdbf;
        [ThreadStatic]
        static double sumdbg;
        [ThreadStatic]
        static double sumdbh;
        [ThreadStatic]
        static double sumdbi;
        [ThreadStatic]
        static double sumdbj;
        [ThreadStatic]
        static double sumdbk;
        [ThreadStatic]
        static double sumdbl;
        [ThreadStatic]
        static double sumdcd;
        [ThreadStatic]
        static double sumdce;
        [ThreadStatic]
        static double sumdcf;
        [ThreadStatic]
        static double sumdcg;
        [ThreadStatic]
        static double sumdch;
        [ThreadStatic]
        static double sumdci;
        [ThreadStatic]
        static double sumdcj;
        [ThreadStatic]
        static double sumdck;
        [ThreadStatic]
        static double sumdcl;
        [ThreadStatic]
        static double sumdde;
        [ThreadStatic]
        static double sumddf;
        [ThreadStatic]
        static double sumddg;
        [ThreadStatic]
        static double sumddh;
        [ThreadStatic]
        static double sumddi;
        [ThreadStatic]
        static double sumddj;
        [ThreadStatic]
        static double sumddk;
        [ThreadStatic]
        static double sumddl;
        [ThreadStatic]
        static double sumdef;
        [ThreadStatic]
        static double sumdeg;
        [ThreadStatic]
        static double sumdeh;
        [ThreadStatic]
        static double sumdei;
        [ThreadStatic]
        static double sumdej;
        [ThreadStatic]
        static double sumdek;
        [ThreadStatic]
        static double sumdel;
        [ThreadStatic]
        static double sumdfg;
        [ThreadStatic]
        static double sumdfh;
        [ThreadStatic]
        static double sumdfi;
        [ThreadStatic]
        static double sumdfj;
        [ThreadStatic]
        static double sumdfk;
        [ThreadStatic]
        static double sumdfl;
        [ThreadStatic]
        static double sumdgh;
        [ThreadStatic]
        static double sumdgi;
        [ThreadStatic]
        static double sumdgj;
        [ThreadStatic]
        static double sumdgk;
        [ThreadStatic]
        static double sumdgl;
        [ThreadStatic]
        static double sumdhi;
        [ThreadStatic]
        static double sumdhj;
        [ThreadStatic]
        static double sumdhk;
        [ThreadStatic]
        static double sumdhl;
        [ThreadStatic]
        static double sumdij;
        [ThreadStatic]
        static double sumdik;
        [ThreadStatic]
        static double sumdil;
        [ThreadStatic]
        static double sumdjk;
        [ThreadStatic]
        static double sumdjl;
        [ThreadStatic]
        static double sumdkl;
        [ThreadStatic]
        static double sumdmb;
        [ThreadStatic]
        static double sumdmc;
        [ThreadStatic]
        static double sumdmd;
        [ThreadStatic]
        static double sumdme;
        [ThreadStatic]
        static double sumdmf;
        [ThreadStatic]
        static double sumdmg;
        [ThreadStatic]
        static double sumdmh;
        [ThreadStatic]
        static double sumdmi;
        [ThreadStatic]
        static double sumdmj;
        [ThreadStatic]
        static double sumdmk;
        [ThreadStatic]
        static double sumdml;
        [ThreadStatic]
        static double sumdma;
        [ThreadStatic]
        static double iagk;
        [ThreadStatic]
        static double iagk2;
        [ThreadStatic]
        static double iagk3;
        [ThreadStatic]
        static double iagk4;
        [ThreadStatic]
        static double iagk5;
        [ThreadStatic]
        static double iagk6;
        [ThreadStatic]
        static double iagk7;
        [ThreadStatic]
        static double iagk8;
        [ThreadStatic]
        static double iagk9;
        [ThreadStatic]
        static double iagk10;
        [ThreadStatic]
        static double iagk11;
        [ThreadStatic]
        static double iagk12;
        [ThreadStatic]
        static double iagk0;
        [ThreadStatic]
        static double[,] H;
        [ThreadStatic]
        static double[,] H1;
        [ThreadStatic]
        static double[,] NR;
        [ThreadStatic]
        static double[,] NRe;
        [ThreadStatic]
        static double[,] SE;
        [ThreadStatic]
        static double[,] V;
        [ThreadStatic]
        static double[,] IV;
        [ThreadStatic]
        static double[,] Z;
        [ThreadStatic]
        static double[,] IH;
        [ThreadStatic]
        static double[,] ini;
        [ThreadStatic]
        static double[,] ini1;
        [ThreadStatic]
        static double[,] ini2;
        [ThreadStatic]
        static double[,] newArray;
        [ThreadStatic]
        static double[,] newArray1;
        [ThreadStatic]
        static double[,] newArray2;
        [ThreadStatic]
        static double[,] newArray3;
        [ThreadStatic]
        static double[,] newArray4;
        [ThreadStatic]
        static double[,] newArray5;
        [ThreadStatic]
        static double z0;
        [ThreadStatic]
        static double z1;
        [ThreadStatic]
        static double z2;
        [ThreadStatic]
        static double z3;
        [ThreadStatic]
        static double z4;
        [ThreadStatic]
        static double z5;
        [ThreadStatic]
        static double z6;
        [ThreadStatic]
        static double z7;
        [ThreadStatic]
        static double z8;
        [ThreadStatic]
        static double z9;
        [ThreadStatic]
        static double z10;
        [ThreadStatic]
        static double z11;
        [ThreadStatic]
        static double z12;
        [ThreadStatic]
        static double[,] zscore;
        [ThreadStatic]
        static double[,] iagkmatrix;
        [ThreadStatic]
        static double[,] waldmatrix;
        [ThreadStatic]
        static double[] SNP;
        [ThreadStatic]
        static double[] Cov1;
        [ThreadStatic]
        static double[] Cov2;
        [ThreadStatic]
        static double[] Cov3;
        [ThreadStatic]
        static double[] Cov4;
        [ThreadStatic]
        static double[] Cov5;
        [ThreadStatic]
        static double[] Cov6;
        [ThreadStatic]
        static double[] Cov7;
        [ThreadStatic]
        static double[] Cov8;
        [ThreadStatic]
        static double[] Cov9;
        [ThreadStatic]
        static double[] Cov10;
        [ThreadStatic]
        static double[] da;
        [ThreadStatic]
        static double[] db;
        [ThreadStatic]
        static double[] dc;
        [ThreadStatic]
        static double[] dd;
        [ThreadStatic]
        static double[] de;
        [ThreadStatic]
        static double[] df;
        [ThreadStatic]
        static double[] dg;
        [ThreadStatic]
        static double[] dh;
        [ThreadStatic]
        static double[] di;
        [ThreadStatic]
        static double[] dj;
        [ThreadStatic]
        static double[] dk;
        [ThreadStatic]
        static double[] dl;
        [ThreadStatic]
        static double[] d2a;
        [ThreadStatic]
        static double[] d2b;
        [ThreadStatic]
        static double[] d2c;
        [ThreadStatic]
        static double[] d2d;
        [ThreadStatic]
        static double[] d2e;
        [ThreadStatic]
        static double[] d2f;
        [ThreadStatic]
        static double[] d2g;
        [ThreadStatic]
        static double[] d2h;
        [ThreadStatic]
        static double[] d2i;
        [ThreadStatic]
        static double[] d2j;
        [ThreadStatic]
        static double[] d2k;
        [ThreadStatic]
        static double[] d2l;
        [ThreadStatic]
        static double[] dab;
        [ThreadStatic]
        static double[] dac;
        [ThreadStatic]
        static double[] dad;
        [ThreadStatic]
        static double[] dae;
        [ThreadStatic]
        static double[] daf;
        [ThreadStatic]
        static double[] dag;
        [ThreadStatic]
        static double[] dah;
        [ThreadStatic]
        static double[] dai;
        [ThreadStatic]
        static double[] daj;
        [ThreadStatic]
        static double[] dak;
        [ThreadStatic]
        static double[] dal;
        [ThreadStatic]
        static double[] dbc;
        [ThreadStatic]
        static double[] dbd;
        [ThreadStatic]
        static double[] dbe;
        [ThreadStatic]
        static double[] dbf;
        [ThreadStatic]
        static double[] dbg;
        [ThreadStatic]
        static double[] dbh;
        [ThreadStatic]
        static double[] dbi;
        [ThreadStatic]
        static double[] dbj;
        [ThreadStatic]
        static double[] dbk;
        [ThreadStatic]
        static double[] dbl;
        [ThreadStatic]
        static double[] dcd;
        [ThreadStatic]
        static double[] dce;
        [ThreadStatic]
        static double[] dcf;
        [ThreadStatic]
        static double[] dcg;
        [ThreadStatic]
        static double[] dch;
        [ThreadStatic]
        static double[] dci;
        [ThreadStatic]
        static double[] dcj;
        [ThreadStatic]
        static double[] dck;
        [ThreadStatic]
        static double[] dcl;
        [ThreadStatic]
        static double[] dde;
        [ThreadStatic]
        static double[] ddf;
        [ThreadStatic]
        static double[] ddg;
        [ThreadStatic]
        static double[] ddh;
        [ThreadStatic]
        static double[] ddi;
        [ThreadStatic]
        static double[] ddj;
        [ThreadStatic]
        static double[] ddk;
        [ThreadStatic]
        static double[] ddl;
        [ThreadStatic]
        static double[] def;
        [ThreadStatic]
        static double[] deg;
        [ThreadStatic]
        static double[] deh;
        [ThreadStatic]
        static double[] dei;
        [ThreadStatic]
        static double[] dej;
        [ThreadStatic]
        static double[] dek;
        [ThreadStatic]
        static double[] del;
        [ThreadStatic]
        static double[] dfg;
        [ThreadStatic]
        static double[] dfh;
        [ThreadStatic]
        static double[] dfi;
        [ThreadStatic]
        static double[] dfj;
        [ThreadStatic]
        static double[] dfk;
        [ThreadStatic]
        static double[] dfl;
        [ThreadStatic]
        static double[] dgh;
        [ThreadStatic]
        static double[] dgi;
        [ThreadStatic]
        static double[] dgj;
        [ThreadStatic]
        static double[] dgk;
        [ThreadStatic]
        static double[] dgl;
        [ThreadStatic]
        static double[] dhi;
        [ThreadStatic]
        static double[] dhj;
        [ThreadStatic]
        static double[] dhk;
        [ThreadStatic]
        static double[] dhl;
        [ThreadStatic]
        static double[] dij;
        [ThreadStatic]
        static double[] dik;
        [ThreadStatic]
        static double[] dil;
        [ThreadStatic]
        static double[] djk;
        [ThreadStatic]
        static double[] djl;
        [ThreadStatic]
        static double[] dkl;
        [ThreadStatic]
        static double[] dm;
        [ThreadStatic]
        static double[] d2m;
        [ThreadStatic]
        static double[] dma;
        [ThreadStatic]
        static double[] dmb;
        [ThreadStatic]
        static double[] dmc;
        [ThreadStatic]
        static double[] dmd;
        [ThreadStatic]
        static double[] dme;
        [ThreadStatic]
        static double[] dmf;
        [ThreadStatic]
        static double[] dmg;
        [ThreadStatic]
        static double[] dmh;
        [ThreadStatic]
        static double[] dmi;
        [ThreadStatic]
        static double[] dmj;
        [ThreadStatic]
        static double[] dmk;
        [ThreadStatic]
        static double[] dml;
        [ThreadStatic]
        static double b0pval;
        [ThreadStatic]
        static double b1pval;
        [ThreadStatic]
        static double b2pval;
        [ThreadStatic]
        static double b3pval;
        [ThreadStatic]
        static double b4pval;
        [ThreadStatic]
        static double b5pval;
        [ThreadStatic]
        static double b6pval;
        [ThreadStatic]
        static double b7pval;
        [ThreadStatic]
        static double b8pval;
        [ThreadStatic]
        static double b9pval;
        [ThreadStatic]
        static double b10pval;
        [ThreadStatic]
        static double b11pval;
        [ThreadStatic]
        static double b12pval;
        [ThreadStatic]
        static WaldTest w0;
        [ThreadStatic]
        static WaldTest w1;
        [ThreadStatic]
        static WaldTest w2;
        [ThreadStatic]
        static WaldTest w3;
        [ThreadStatic]
        static WaldTest w4;
        [ThreadStatic]
        static WaldTest w5;
        [ThreadStatic]
        static WaldTest w6;
        [ThreadStatic]
        static WaldTest w7;
        [ThreadStatic]
        static WaldTest w8;
        [ThreadStatic]
        static WaldTest w9;
        [ThreadStatic]
        static WaldTest w10;
        [ThreadStatic]
        static WaldTest w11;
        [ThreadStatic]
        static WaldTest w12;
        [ThreadStatic]
        static int covcount;
        [ThreadStatic]
        static double loglikenull;
        [ThreadStatic]
        static double loglikealt;
        [ThreadStatic]
        static double likeratiotest;
        [ThreadStatic]
        static ChiSquareTest chi;
        [ThreadStatic]
        static double jointint;
        [ThreadStatic]
        static double[] weiblike;
        [ThreadStatic]
        static double[] weiblike2;
        private static SemaphoreSlim _semaSlim;

        static void Main(string[] args)
        {

            Console.WriteLine("");
            Console.WriteLine("                              Hello, Welcome to SurvivalGWAS_SV");
            Console.WriteLine("$-------------------------------------------------------------------------------------------$");
            Console.WriteLine("|                                       June, 2016                                          |");
            Console.WriteLine("|-------------------------------------------------------------------------------------------|");
            Console.WriteLine("|                  (C) 2016 Hamzah Syed, Andrea L Jorgensen & Andrew P Morris               |");
            Console.WriteLine("|                               GNU General Public License, v3                              |");
            Console.WriteLine("|-------------------------------------------------------------------------------------------|");
            Console.WriteLine("|       SurvivalGWAS_SV - Genome-wide association study analysis of imputed genotypes       |");
            Console.WriteLine("|                                 with time-to-event outcomes                               |");
            Console.WriteLine("|                                                                                           |");
            Console.WriteLine("|             This single variant analytics tool is part of the SurvivalGWAS Suite          |");
            Console.WriteLine("|                                                                                           |");
            Console.WriteLine("|                 For documentation, citation & bug-report instructions:                    |");
            Console.WriteLine("|https://www.liverpool.ac.uk/translational-medicine/research/statistical-genetics/software/ |");
            Console.WriteLine("$-------------------------------------------------------------------------------------------$");
            Console.WriteLine("");

            // Commands
            var p = new OptionSet() {
            {"gf|gen_file=", "The name of the genotype file",
              (string v) => file = v },
            { "sf|sample_file=", "The name of the sample file.",
              (string v) => file2 = v },
            { "t|time=", "The observation time",
              (string v) => Stime = v },
            { "c|censor=", "The censoring indicator",
              (string v) => cens = v },
            { "cov|covariates=", "A list of covariates. Each one seperated by a comma (,)",
              (string v) => covariates = v },
            { "i|int=", "The interaction between SNP and one covariate. Seperate using a comma (,)",
              (string v) => int1 = v },
            { "lstart|linestart=", "Specify line in file start position for more efficient program runtime",
              (string v) => num = v },
            { "lstop|linestop=", "Specify line in file stop position for more efficient program runtime",
              (string v) => num2 = v },
            { "sp|start_position=", "The start position on the chromosome. You still need to specify the number of lines in file using -lstart & -lstop commands",
              (string v) => start = v },
            { "ep|end_position=", "The stop position on the chromosome. -sp & -ep commands are substantially slower than using the -lstart & -lstop on their own",
              (string v) => stop = v },
            { "chr|chromosome=", "User specified chromosome number",
              (string v) => chrome = v },
            { "p|print=", "Enter 'onlysnp' if you want only the SNP analysis output to be in the output file and 'onlyint' if you want only the interaction analysis output to be in the output file",
              (string v) => print = v },
            { "m|method=", "Specify choice of method for analysis",
              (string v) => method = v },
            { "o|output=", "Name of file for output to be saved in",
              (string v) => SaveFile = v },
            { "threads|cpu=", "Number of threads. On a multi-core system, multiple threads can execute tasks in parallel, with each core executing a different thread",
              (string v) => cores = v },
            { "h|help",  "Command Help",
              v => show_help = v != null },
            };

            List<string> extra;
            try
            {
                extra = p.Parse(args);
            }
            catch (OptionException e)
            {
                Console.WriteLine(e.Message);
                Console.WriteLine("Try '-help' for more information.");
                return;
            }

            if (show_help)
            {
                ShowHelp(p);
                return;
            }

            // Specify directory and file names, with validation

            if (File.Exists(file) == false || File.Exists(file2) == false)
            {
                Console.WriteLine("This input file does not exist.");
                System.Environment.Exit(1);
            }

            // Validate threads input

            if (cores == null || cores == "0")
            {
                Console.WriteLine("Please specify number of threads (-threads=) to be created");
                System.Environment.Exit(1);
            }

            // Validate method input

            if (method == null || (method != "cox" && method != "weibull"))
            {
                Console.WriteLine("Incorrect method (-m=) entered");
                System.Environment.Exit(1);
            }
            // Validate print option name

            if (print != null && (print != "onlysnp" && print != "onlyint"))
            {
                Console.WriteLine("The input for the -print (-p=) option is incorrect.");
                System.Environment.Exit(1);
            }
            // User defines range by lines in file using num and num2 or chromosome position using start and stop 
            if (num == null || num2 == null) // || (start != null && stop != null && num != null && num2 != null))
            {
                Console.WriteLine("Number of lines in file or chromosome position interval to analyse not specified. Please use -help command.");
                System.Environment.Exit(1);
            }

            // Create output file with heading

            if (SaveFile == null)
            {
                Console.WriteLine("File to save output not specified. Use -o= command");
                System.Environment.Exit(1);
            }
            else
            {
                if (method == "cox")
                {
                    if (File.Exists(SaveFile))
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                        }
                        Console.WriteLine(SaveFile + " already exists in directory, new lines have been added to existing file");
                    }
                    else
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            writer.WriteLine("InputName rsid Chr Pos EA NonEA CoefValue HR SE LowerCI UpperCI Waldpv LRTpv ModLRTpv Jointassoc EAF MAF Infoscore");
                        }
                    }
                }
                else
                {
                    if (File.Exists(SaveFile))
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                        }
                        Console.WriteLine(SaveFile + " already exists in directory, new lines have been added to existing file");
                    }
                    else
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            writer.WriteLine("InputName rsid Chr Pos EA NonEA CoefValue AF HR(PH) SE zscorestat p-value Waldpv Jointassoc EAF MAF Infoscore Shape");
                        }
                    }
                }
            }
            _semaSlim = new SemaphoreSlim(int.Parse(cores)); //user specified number
            Console.WriteLine("Reading genotype data....");

            indexes = new List<long>();

            // File read in (zipped or unzipped)

            if (file.EndsWith(".zip"))
            {
                using (ZipArchive archive = ZipFile.OpenRead(file))
                {
                    foreach (ZipArchiveEntry entry in archive.Entries)
                    {
                        try
                        {
                            if (entry.FullName.EndsWith(".gen", StringComparison.OrdinalIgnoreCase))
                            {
                                entry.ExtractToFile(Path.Combine(Directory.GetCurrentDirectory(), entry.FullName));
                            }
                            else if (entry.FullName.EndsWith(".impute", StringComparison.OrdinalIgnoreCase))
                            {
                                entry.ExtractToFile(Path.Combine(Directory.GetCurrentDirectory(), entry.FullName));
                            }
                            else if (entry.FullName.EndsWith(".vcf", StringComparison.OrdinalIgnoreCase))
                            {
                                entry.ExtractToFile(Path.Combine(Directory.GetCurrentDirectory(), entry.FullName));
                            }
                        }
                        catch (Exception d)
                        {
                            Console.WriteLine("File has already been unzipped in this directory");
                        }

                        var fs = File.OpenRead(entry.FullName);
                        indexes.Add(fs.Position);

                        while ((chr = fs.ReadByte()) != -1)
                        {
                            if (chr == '\n')
                            {
                                indexes.Add(fs.Position);
                            }
                        }
                    }
                }
            }
            else if (file.EndsWith(".gz"))
            {
                FileInfo gzipFileName = new FileInfo(file);
                using (FileStream fileToDecompressAsStream = gzipFileName.OpenRead())
                {
                    decompressedFileName = file.Remove(file.Length - 3);
                    using (FileStream decompressedStream = File.Create(decompressedFileName))
                    {
                        using (GZipStream decompressionStream = new GZipStream(fileToDecompressAsStream, CompressionMode.Decompress))
                        {
                            try
                            {
                                decompressionStream.CopyTo(decompressedStream);
                            }
                            catch (Exception ex)
                            {
                            }
                        }
                    }
                    var fs = File.OpenRead(decompressedFileName);
                    indexes.Add(fs.Position);
                    while ((chr = fs.ReadByte()) != -1)
                    {
                        if (chr == '\n')
                        {
                            indexes.Add(fs.Position);
                        }
                    }
                }
            }
            else
            {
                using (var fs = File.OpenRead(file))
                {
                    indexes.Add(fs.Position);
                    while ((chr = fs.ReadByte()) != -1)
                    {
                        if (chr == '\n')
                        {
                            indexes.Add(fs.Position);
                        }
                    }
                }
            }

            Console.WriteLine("Data read in complete");

            // initialising threads  
            for (int i = 1; i <= int.Parse(cores); i++)
            {
                new Thread(DoWork).Start(i);
            }

        }

        public static void DoWork(object threadId)
        {
            _semaSlim.Wait();
            Console.WriteLine("Batch {0} has started", threadId);

            int lineCount = indexes.Count;
            string lineContent = "";
            // num2 == "max" finds total number of lines in file. This helps with for loops in shell scripts.
            if (threadId.Equals(1))
            {
                number = int.Parse(num);
                if (num2 == "max")
                {
                    number2 = ((indexes.Count - int.Parse(num)) / int.Parse(cores)) + number;
                }
                else
                {
                    number2 = ((int.Parse(num2) - int.Parse(num)) / int.Parse(cores)) + number;
                }
            }
            else
            {
                if (num2 == "max")
                {
                    number = (((indexes.Count - int.Parse(num)) / int.Parse(cores)) * (Convert.ToInt32(threadId) - 1)) + 1 + int.Parse(num);
                    number2 = (((indexes.Count - int.Parse(num)) / int.Parse(cores)) * Convert.ToInt32(threadId)) + int.Parse(num);
                }
                else
                {

                    number = (((int.Parse(num2) - int.Parse(num)) / int.Parse(cores)) * (Convert.ToInt32(threadId) - 1)) + 1 + int.Parse(num);
                    number2 = (((int.Parse(num2) - int.Parse(num)) / int.Parse(cores)) * Convert.ToInt32(threadId)) + int.Parse(num);
                }
            }
          
            int randLineNum = number;
            int counter = 0;

            while (lineContent != null)
            {
                // Read the line in genotype file
                if (file.EndsWith(".zip"))
                {
                    using (ZipArchive archive = ZipFile.OpenRead(file))
                    {
                        foreach (ZipArchiveEntry entry2 in archive.Entries)
                        {
                            using (var fs = File.OpenRead(entry2.FullName))
                            {
                                fs.Position = indexes[randLineNum];
                                using (var sr = new StreamReader(fs))
                                {
                                    lineContent = sr.ReadLine();
                                }
                            }
                        }
                        if (lineContent == null)
                        {
                            break;
                        }
                    }
                }
                else if (file.EndsWith(".gz"))
                {
                    using (var fs = File.OpenRead(decompressedFileName))
                    {
                        fs.Position = indexes[randLineNum];
                        using (var sr = new StreamReader(fs))
                        {
                            lineContent = sr.ReadLine();
                        }
                    }
                    if (lineContent == null)
                    {
                        break;
                    }
                }
                else
                {
                    using (var fs = File.OpenRead(file))
                    {
                        fs.Position = indexes[randLineNum];
                        using (var sr = new StreamReader(fs))
                        {
                            lineContent = sr.ReadLine();
                        }
                    }
                    if (lineContent == null)
                    {
                        break;
                    }
                }

                table = new DataTable();

                //reads gen file data into table
                if (file.Contains(".gen") || file.Contains(".impute"))
                {
                    if (lineContent == "")
                    {
                        goto FINISHED;
                    }
                    splitFileContents2 = lineContent.Split(' ').ToArray();
                    a = splitFileContents2.GetRows(0)[0];
                    table.Columns.Add(a[0]);

                    for (int i = 2; i < splitFileContents2.Count(); i++)
                    {
                        DataRow row = table.NewRow();
                        row.ItemArray = splitFileContents2[i];
                        table.Rows.Add(row);
                    }
                    // Stopping criteria for gen file

                    if ((double.Parse(table.Rows[0][0].ToString()) > Convert.ToDouble(stop) && randLineNum > number2) || randLineNum > number2)
                    {
                        goto FINISHED;
                    }

                    if ((double.Parse(table.Rows[0][0].ToString()) >= Convert.ToDouble(start) && double.Parse(table.Rows[0][0].ToString()) <= Convert.ToDouble(stop) && randLineNum >= number && randLineNum <= number2) || (randLineNum >= number && randLineNum <= number2 && start == null & stop == null))
                    {

                        // gen file conversion, probabilities.                       
                        dos = new DataTable();
                        NN = new List<double>();
                        PP = new List<double>();
                        FF = new List<double>();
                        FPsub = new List<double>();

                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();

                        dos = table.Clone();

                        for (int i = 0; i < (table.Rows.Count) / 3; i++)
                        {
                            // convert probabilities to dosages using additive model

                            DataRow row2 = dos.NewRow();
                            dos.Rows.Add(row2);

                            dos.Rows[i][0] = Convert.ToDouble(table.Rows[3 * i + 1][0]) + 2 * (Convert.ToDouble(table.Rows[3 * i + 2][0]));

                            // allele frequency calculation

                            double N = Convert.ToDouble(table.Rows[3 * i][0]) + Convert.ToDouble(table.Rows[3 * i + 1][0]) + Convert.ToDouble(table.Rows[3 * i + 2][0]);
                            NN.Add(N);

                            double P = Convert.ToDouble(table.Rows[3 * i + 1][0]) + 2 * (Convert.ToDouble(table.Rows[3 * i + 2][0]));
                            PP.Add(P);

                            double F = Convert.ToDouble(table.Rows[3 * i + 1][0]) + 4 * (Convert.ToDouble(table.Rows[3 * i + 2][0]));
                            FF.Add(F);

                            double FP = F - Math.Pow(P, 2);
                            FPsub.Add(FP);
                        }

                        double sumNN = NN.Sum();
                        double sumPP = PP.Sum();
                        double sumFP = FPsub.Sum();
                        allele = sumPP / (2 * sumNN); // population allele frequency
                        infoscore = (((2 * sumNN) / (allele * (1 - allele))) - ((sumFP) / (Math.Pow(allele, 2) * Math.Pow(1 - allele, 2)))) / ((2 * sumNN) / (allele * (1 - allele)));

                        if (allele > 0.5)
                        {
                            allele = 1 - allele; // MAF                               
                        }
                    }
                    else
                    {
                        goto CRITFAIL;
                    }
                }
                else
                {
                    //reads vcf file data in table

                    if (lineContent == "")
                    {
                        goto FINISHED;
                    }
                    else if (lineContent.StartsWith("#") || Char.IsLetter(lineContent[0]))
                    {
                        goto CRITFAIL;
                    }

                    splitFileContents2 = lineContent.Split('\t', ' ').ToArray();

                    // this takes the id which is 3rd item
                    a = splitFileContents2.GetRows(2)[0];
                    table.Columns.Add(a[0]);

                    for (int i = 0; i < splitFileContents2.Count(); i++)
                    {
                        DataRow row = table.NewRow();
                        row.ItemArray = splitFileContents2[i];
                        table.Rows.Add(row);
                    }

                    if ((double.Parse(table.Rows[1][0].ToString()) > Convert.ToDouble(stop) && randLineNum > number2) || randLineNum > number2)
                    {
                        goto FINISHED;
                    }

                    if (((double.Parse(table.Rows[1][0].ToString()) >= Convert.ToDouble(start)) && (double.Parse(table.Rows[1][0].ToString()) <= Convert.ToDouble(stop)) && (randLineNum >= number && randLineNum <= number2)) || (randLineNum >= number && randLineNum <= number2 && start == null & stop == null))
                    {

                        dos = new DataTable();
                        NN = new List<double>();
                        PP = new List<double>();
                        FF = new List<double>();
                        FPsub = new List<double>();
                        gsum2 = new List<double>();

                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();
                        table.Rows.RemoveAt(0);
                        table.AcceptChanges();

                        dos = table.Clone();

                        // Hard genotype if dosages or probablities not available.
                        // convert vcf 0|0, 0|1, 1|1 to dosages
                        if (!splitFileContents2[8][0].Contains("GP") && !splitFileContents2[8][0].Contains("DS") && splitFileContents2[8][0].Contains("GT"))
                        {
                            for (int i = 0; i < (table.Rows.Count); i++)
                            {
                                DataRow row2 = dos.NewRow();
                                dos.Rows.Add(row2);

                                double genotype = 0.0;

                                if (table.Rows[i][0].ToString().StartsWith("0|0") || table.Rows[i][0].ToString().StartsWith("0/0"))
                                {
                                    genotype = 0;
                                }
                                else if (table.Rows[i][0].ToString().StartsWith("0|1") || table.Rows[i][0].ToString().StartsWith("0/1")) //can be 1/0
                                {
                                    genotype = 1;
                                }
                                else if (table.Rows[i][0].ToString().StartsWith("1|1") || table.Rows[i][0].ToString().StartsWith("1/1"))
                                {
                                    genotype = 2;
                                }

                                dos.Rows[i][0] = genotype;

                                double gsum = Convert.ToInt32(dos.Rows[i][0]);
                                gsum2.Add(gsum);
                            }

                            // extract maf from file
                            string[] vcfsplit = splitFileContents2[7][0].Split(';');
                            int mafnum;
                            if (vcfsplit.Contains("MAF"))
                            {
                                mafnum = Array.IndexOf(vcfsplit, "MAF");
                                string[] vcfsplit2 = vcfsplit[mafnum].Split('=');
                                double vcfmaf = Convert.ToDouble(vcfsplit2[1]);
                                allele = vcfmaf;
                            }
                            else if (vcfsplit.Contains("AF"))
                            {
                                mafnum = Array.IndexOf(vcfsplit, "AF");
                                string[] vcfsplit2 = vcfsplit[mafnum].Split('=');
                                double vcfmaf = Convert.ToDouble(vcfsplit2[1]);
                                allele = vcfmaf;
                            }
                            else
                            {

                                double sumgeno = gsum2.Sum();

                                allele = sumgeno / (2 * dos.Rows.Count);
                            }

                            if (allele > 0.5)
                            {
                                allele = 1 - allele; // MAF
                            }
                            infoscore = -1;
                        }

                        else if (splitFileContents2[8][0].Contains("DS") && !splitFileContents2[8][0].Contains("GP"))
                        {
                            for (int i = 0; i < (table.Rows.Count); i++)
                            {

                                DataRow row2 = dos.NewRow();
                                dos.Rows.Add(row2);

                                string[] ds = table.Rows[i][0].ToString().Split(':');
                                string[] dsrow = splitFileContents2[8][0].Split(':');
                                int dsnum = Array.IndexOf(dsrow, "DS");
                                string ds2 = ds[dsnum];
                                dos.Rows[i][0] = Convert.ToDouble(ds2);
                                double P = Convert.ToDouble(ds2);
                                PP.Add(P);
                            }
                            double sumNN = dos.Rows.Count;
                            double sumPP = PP.Sum();
                            allele = sumPP / (2 * sumNN); // population allele frequency
                            if (allele > 0.5)
                            {
                                allele = 1 - allele; // MAF
                            }
                            infoscore = -1;
                        }

                        else if (splitFileContents2[8][0].Contains("GP"))
                        {
                            for (int i = 0; i < (table.Rows.Count); i++)
                            {

                                DataRow row2 = dos.NewRow();
                                dos.Rows.Add(row2);

                                string[] gp = table.Rows[i][0].ToString().Split(':');
                                string[] gprow = splitFileContents2[8][0].Split(':');
                                int gpnum = Array.IndexOf(gprow, "GP");
                                string[] gp2 = gp[gpnum].Split(',');

                                dos.Rows[i][0] = Convert.ToDouble(gp2[1]) + 2 * (Convert.ToDouble(gp2[2]));

                                double N = Convert.ToDouble(gp2[0]) + Convert.ToDouble(gp2[1]) + Convert.ToDouble(gp2[2]);
                                NN.Add(N);

                                double P = Convert.ToDouble(gp2[1]) + 2 * (Convert.ToDouble(gp2[2]));
                                PP.Add(P);

                                double F = Convert.ToDouble(gp2[1]) + 4 * (Convert.ToDouble(gp2[2]));
                                FF.Add(F);

                                double FP = F - Math.Pow(P, 2);
                                FPsub.Add(FP);
                            }
                            // this point onwards it is same as above
                            // allele frequency calculation
                            double sumNN = NN.Sum();
                            double sumPP = PP.Sum();
                            double sumFP = FPsub.Sum();
                            allele = sumPP / (2 * sumNN); // population allele frequency
                            infoscore = (((2 * sumNN) / (allele * (1 - allele))) - ((sumFP) / (Math.Pow(allele, 2) * Math.Pow(1 - allele, 2)))) / ((2 * sumNN) / (allele * (1 - allele)));

                            if (allele > 0.5)
                            {
                                allele = 1 - allele; // MAF                               
                            }
                        }
                    }
                    else
                    {
                        goto CRITFAIL;
                    }
                }


                // sample file read in
                dt = new DataTable();
                independent = new DataTable();
                names = new List<string>();
                
                string[] columns = null;

                var lines = File.ReadAllLines(file2);

                // assuming the first row contains the columns information
                if (lines.Count() > 0)
                {
                    columns = lines[0].Split(new char[] { ' ', '\t' });

                    foreach (var column in columns)
                        dt.Columns.Add(column);
                }

                // reading rest of the data
                for (int i = 1; i < lines.Count(); i++)
                {
                    dr = dt.NewRow();
                    values = lines[i].Split(new char[] { ' ', '\t' });

                    // Converting NA values to NaN
              /*      if (values.Contains("NA"))
                    {
                        for (int z = 0; z < values.Count(); z++)
                        {
                            if (values[z].Contains("NA"))
                            {
                                values[z] = Convert.ToString(Double.NaN);
                            }
                        }
                    }*/

                    for (int j = 0; j < values.Count() && j < columns.Count(); j++)
                        dr[j] = values[j];

                    dt.Rows.Add(dr);
                }

                if (dt.Rows[0][0].Equals("0"))
                {
                    dt.Rows[0].Delete();
                    dt.AcceptChanges();
                }

                //covariates

                if (covariates != null)
                {
                    split = covariates.Split(new Char[] { ',' });

                    for (int g = 0; g < split.Length; g++)
                    {
                        if (dt.Columns.Contains(split[g]))
                        {
                            continue;

                        }
                        else
                        {
                            Console.WriteLine("Column heading " + split[g] + " does not exist!");
                            System.Environment.Exit(1);
                        }
                    }
                }

                //interaction between SNP and covariate

                if (int1 != null)
                {
                    split2 = int1.Split(new Char[] { ',' });

                    for (int v = 0; v < split2.Length; v++)
                    {
                        if (dt.Columns.Contains(split2[1]))
                        {
                            continue;

                        }

                        else
                        {
                            Console.WriteLine("Column heading " + split2[1] + " does not exist!");
                            System.Environment.Exit(1);
                        }
                    }
                }


                if (covariates != null && int1 == null)
                {
                    for (int f = 0; f < dt.Rows.Count; f++)
                    {
                        if ((dt.Rows[f][Stime].Equals("NA")) || (dt.Rows[f][Stime].Equals(0.ToString())))
                        {//Double.NaN.ToString()
                            dt.Rows[f].Delete();
                            dt.AcceptChanges();
                            dos.Rows[f].Delete();
                            dos.AcceptChanges();
                            f = f - 1;
                        }
                    }

                    foreach (string name in split)
                    {

                        names.Add(name);

                        // Deleting rows with missing data
                        for (int e = 0; e < dt.Rows.Count; e++)
                        {
                            if (dt.Rows[e][name].Equals("NA"))
                            {
                                dt.Rows[e].Delete();
                                dt.AcceptChanges();
                                dos.Rows[e].Delete();
                                dos.AcceptChanges();
                                e = e - 1;
                            }
                        }
                       
                    }
                    foreach (string name in split)
                    {
                        independent.Columns.Add(dt.Columns[name].ColumnName);

                        for (int i = 0; i < dos.Rows.Count; i++)
                        {
                            if (name == split[0])
                            {
                                independent.Rows.Add();
                                independent.Rows[i][dt.Columns[name].ColumnName] = dt.Rows[i][name];
                            }
                            else
                            {
                                independent.Rows[i][dt.Columns[name].ColumnName] = dt.Rows[i][name];
                            }
                        }
                    }
                }
                else if (covariates != null && int1 != null)
                {
                    // Deleting rows with missing data
                    for (int f = 0; f < dt.Rows.Count; f++)
                    {
                        if ((dt.Rows[f][Stime].Equals("NA")) || (dt.Rows[f][Stime].Equals(0.ToString())))
                        {
                            dt.Rows[f].Delete();
                            dt.AcceptChanges();
                            dos.Rows[f].Delete();
                            dos.AcceptChanges();
                            f = f - 1;
                        }
                    }

                    foreach (string name in split)
                    {

                        names.Add(name);

                        for (int e = 0; e < dt.Rows.Count; e++)
                        {
                            if (dt.Rows[e][name].Equals("NA"))
                            {
                                dt.Rows[e].Delete();
                                dt.AcceptChanges();
                                dos.Rows[e].Delete();
                                dos.AcceptChanges();
                                e = e - 1;

                            }
                        }                     
                    }
                    foreach (string name in split)
                    {
                        independent.Columns.Add(dt.Columns[name].ColumnName);

                        for (int i = 0; i < dos.Rows.Count; i++)
                        {
                            if (name == split[0])
                            {
                                independent.Rows.Add();
                                independent.Rows[i][dt.Columns[name].ColumnName] = dt.Rows[i][name];
                            }
                            else
                            {
                                independent.Rows[i][dt.Columns[name].ColumnName] = dt.Rows[i][name];
                            }
                        }
                    }
                    intColumn = new DataColumn(dos.Columns[0] + split2[1]);
                    dt.Columns.Add(intColumn);
                    for (int h = 0; h < dt.Rows.Count; h++)
                    {
                        dt.Rows[h][intColumn] = double.Parse(dos.Rows[h][0].ToString()) * double.Parse(dt.Rows[h][split2[1]].ToString());
                    }
                    names.Add(intColumn.ColumnName);

                    independent.Columns.Add(dt.Columns[intColumn.ColumnName].ColumnName);

                    for (int i = 0; i < dos.Rows.Count; i++)
                    {
                        independent.Rows[i][dt.Columns[intColumn.ColumnName].ColumnName] = dt.Rows[i][intColumn];
                    }
                }

                else
                {
                    for (int f = 0; f < dt.Rows.Count; f++)
                    {
                        if ((dt.Rows[f][Stime].Equals("NA")) || (dt.Rows[f][Stime].Equals(0.ToString())))
                        {
                            dt.Rows[f].Delete();
                            dt.AcceptChanges();
                            dos.Rows[f].Delete();
                            dos.AcceptChanges();
                            f = f - 1;
                        }
                    }
                }

                names.Add(dos.Columns[0].ColumnName);
                newnames = names.ToArray();

                independent.Columns.Add(dos.Columns[0].ColumnName);

                for (int i = 0; i < dos.Rows.Count; i++)
                {
                    if (covariates == null)
                    {
                        independent.Rows.Add();
                    }
                    independent.Rows[i][dos.Columns[0].ColumnName] = dos.Rows[i][0];
                }

                //time variable

                if (dt.Columns.Contains(Stime))
                {
                    time = dt.Columns[Stime].ToArray();
                }
                else
                {
                    Console.WriteLine("Column heading " + Stime + " does not exist!");
                    System.Environment.Exit(1);
                }


                //censoring indicator

                if (dt.Columns.Contains(cens))
                {
                    censor = dt.Columns[cens].ToArray().ToInt32();
                }
                else
                {
                    Console.WriteLine("Column heading " + cens + " does not exist!");
                    System.Environment.Exit(1);
                }

                //input variables
                inputer = independent.ToArray();

                if (method == "cox")
                {
                    //cox proportional hazards model

                    pha = new ProportionalHazardsAnalysis(inputer, time, censor, newnames, "", "");

                    try
                    {
                        pha.Compute();
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine("NaN's produced, " + dos.Columns[0] + " analysis did not converge");
                        goto SKIP;
                    }

                    //only if int is selected
                    //null h(t)=h_0(t)exp(b2*cov+...bn*covn)
                    if (covariates != null & int1 != null)
                    {
                        covcount = newnames.Count();
                        independent.Columns.RemoveAt(newnames.Count() - 1);
                        independent.AcceptChanges();
                        independent.Columns.RemoveAt(newnames.Count() - 2);
                        independent.AcceptChanges();
                        inputer = independent.ToArray();                        
                        newnames.RemoveAt(names.Count()-1);
                        newnames.RemoveAt(names.Count()-2);

                        pha2 = new ProportionalHazardsAnalysis(inputer, time, censor, newnames, "", "");
                        try
                        {
                            pha2.Compute();

                        }
                        catch (Exception ex)
                        {
                        }
                        loglikenull = pha2.LogLikelihood;
                        loglikealt = pha.LogLikelihood;
                        likeratiotest = 2 * (loglikealt - loglikenull);
                        chi = new ChiSquareTest(likeratiotest,2);
                        jointint = chi.PValue;
                    }
                    rwl.AcquireWriterLock(Timeout.Infinite);
                    //     lock (mylock)
                    //     {
                    if (print == "onlysnp" && newnames.Count() > 1)
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            for (int b = newnames.Count() - 1; b < newnames.Count(); b++)
                            {
                                lrt = "NA";
                                /*  if (newnames.Count() < 3)
                                  {
                                      lrt = "NA";
                                  }
                                  else
                                  {
                                      lrt = Convert.ToString(pha.LikelihoodRatioTests[b]);
                                  }*/
                                if (int1==null)
                                {
                                    jointint = -1;
                                }
                                if (chrome==null & file.Contains(".vcf"))
                                {
                                    chrome = splitFileContents2[0][0];
                                }
                                writer.WriteLine(pha.Coefficients[b].Name + " " + splitFileContents2[1][0] + " " + chrome + " " + splitFileContents2[2][0] + " " + splitFileContents2[3][0] + " " + splitFileContents2[4][0] + " " + Math.Round(pha.Coefficients[b].Value, 8) + " " + Math.Round(pha.Coefficients[b].HazardRatio, 8) + " " + Math.Round(pha.Coefficients[b].StandardError, 8) + " " + Math.Round(pha.Coefficients[b].ConfidenceLower, 8) + " " + Math.Round(pha.Coefficients[b].ConfidenceUpper, 8) + " " + pha.Coefficients[b].Wald.PValue + " " + lrt + " " + pha.ChiSquare.PValue + " " + jointint + " " + Math.Round(1 - allele, 8).ToString() + " " + Math.Round(allele, 8).ToString() + " " + Math.Round(infoscore, 8).ToString());
                            }

                        }

                    }
                    else if (print == "onlyint" && newnames.Count() > 1)
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            for (int b = newnames.Count() - 2; b < newnames.Count() - 1; b++)
                            {
                                lrt = "NA";
                                /*  if (newnames.Count() < 3)
                                  {
                                      lrt = "NA";
                                  }
                                  else
                                  {
                                      lrt = Convert.ToString(pha.LikelihoodRatioTests[b]);
                                  }*/
                                if (chrome == null & file.Contains(".vcf"))
                                {
                                    chrome = splitFileContents2[0][0];
                                }
                                writer.WriteLine(pha.Coefficients[b].Name + " " + splitFileContents2[1][0] + " " + chrome + " " + splitFileContents2[2][0] + " " + splitFileContents2[3][0] + " " + splitFileContents2[4][0] + " " + Math.Round(pha.Coefficients[b].Value, 8) + " " + Math.Round(pha.Coefficients[b].HazardRatio, 8) + " " + Math.Round(pha.Coefficients[b].StandardError, 8) + " " + Math.Round(pha.Coefficients[b].ConfidenceLower, 8) + " " + Math.Round(pha.Coefficients[b].ConfidenceUpper, 8) + " " + pha.Coefficients[b].Wald.PValue + " " + lrt + " " + pha.ChiSquare.PValue + " " + jointint + " " + Math.Round(1 - allele, 8).ToString() + " " + Math.Round(allele, 8).ToString() + " " + Math.Round(infoscore, 8).ToString());
                            }

                        }
                    }
                    else
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            for (int b = 0; b < newnames.Count(); b++)
                            {
                                lrt = "NA";
                                /*  if (newnames.Count() < 3)
                                  {
                                      lrt = "NA";
                                  }
                                  else
                                  {
                                      lrt = Convert.ToString(pha.LikelihoodRatioTests[b]);
                                  }*/
                                if (int1 == null)
                                {
                                    jointint = -1;
                                }
                                if (chrome == null & file.Contains(".vcf"))
                                {
                                    chrome = splitFileContents2[0][0];
                                }
                                writer.WriteLine(pha.Coefficients[b].Name + " " + splitFileContents2[1][0] + " " + chrome + " " + splitFileContents2[2][0] + " " + splitFileContents2[3][0] + " " + splitFileContents2[4][0] + " " + Math.Round(pha.Coefficients[b].Value, 8) + " " + Math.Round(pha.Coefficients[b].HazardRatio, 8) + " " + Math.Round(pha.Coefficients[b].StandardError, 8) + " " + Math.Round(pha.Coefficients[b].ConfidenceLower, 8) + " " + Math.Round(pha.Coefficients[b].ConfidenceUpper, 8) + " " + pha.Coefficients[b].Wald.PValue + " " + lrt + " " + pha.ChiSquare.PValue + " " + jointint + " " + Math.Round(1 - allele, 8).ToString() + " " + Math.Round(allele, 8).ToString() + " " + Math.Round(infoscore, 8).ToString());
                            }

                        }
                    }
                    //     }
                    rwl.ReleaseWriterLock();
                }

                if (method == "weibull")
                {
                    // weibull names
                    names2 = new List<string>();
                    names2.Add("Intercept");
                    names2.Add(dos.Columns[0].ColumnName);
                    if (covariates != null)
                    {
                        foreach (string name in split)
                        {

                            names2.Add(name);
                        }

                        try
                        {
                            names2.Add(intColumn.ColumnName);
                        }
                        catch (Exception intcol)
                        {
                        }
                    }
                    names2.Add("Shape");
                    SNP = dos.Columns[0].ToArray();  //.ToInt32();
                    Cov1 = new double[dos.Rows.Count];
                    Cov2 = new double[dos.Rows.Count];
                    Cov3 = new double[dos.Rows.Count];
                    Cov4 = new double[dos.Rows.Count];
                    Cov5 = new double[dos.Rows.Count];
                    Cov6 = new double[dos.Rows.Count];
                    Cov7 = new double[dos.Rows.Count];
                    Cov8 = new double[dos.Rows.Count];
                    Cov9 = new double[dos.Rows.Count];
                    Cov10 = new double[dos.Rows.Count];

                    if (covariates != null)
                    {
                        Cov1 = dt.Columns[split[0]].ToArray();

                        if (covariates.Contains(","))
                        {
                            Cov2 = dt.Columns[split[1]].ToArray();
                        }
                        else if (int1 != null && !covariates.Contains(","))
                        {
                            Cov2 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }

                        if (covariates.Count(c => c == ',') >= 2)
                        {
                            Cov3 = dt.Columns[split[2]].ToArray();
                        }

                        else if (int1 != null && covariates.Count(c => c == ',') == 1)
                        {
                            Cov3 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }

                        if (covariates.Count(c => c == ',') >= 3)
                        {
                            Cov4 = dt.Columns[split[3]].ToArray();
                        }
                        else if (int1 != null && covariates.Count(c => c == ',') == 2)
                        {
                            Cov4 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }

                        if (covariates.Count(c => c == ',') >= 4)
                        {
                            Cov5 = dt.Columns[split[4]].ToArray();

                        }
                        else if (int1 != null && covariates.Count(c => c == ',') == 3)
                        {
                            Cov5 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }
                        if (covariates.Count(c => c == ',') >= 5)
                        {
                            Cov6 = dt.Columns[split[5]].ToArray();

                        }
                        else if (int1 != null && covariates.Count(c => c == ',') == 4)
                        {
                            Cov6 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }
                        if (covariates.Count(c => c == ',') >= 6)
                        {
                            Cov7 = dt.Columns[split[6]].ToArray();

                        }
                        else if (int1 != null && covariates.Count(c => c == ',') == 5)
                        {
                            Cov7 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }
                        if (covariates.Count(c => c == ',') >= 7)
                        {
                            Cov8 = dt.Columns[split[7]].ToArray();

                        }
                        else if (int1 != null && covariates.Count(c => c == ',') == 6)
                        {
                            Cov8 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }
                        if (covariates.Count(c => c == ',') >= 8)
                        {
                            Cov9 = dt.Columns[split[8]].ToArray();

                        }
                        else if (int1 != null && covariates.Count(c => c == ',') == 7)
                        {
                            Cov9 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }
                        if (covariates.Count(c => c == ',') >= 9)
                        {
                            Cov10 = dt.Columns[split[9]].ToArray();

                        }
                        else if (int1 != null && covariates.Count(c => c == ',') == 8)
                        {
                            Cov10 = dt.Columns[dos.Columns[0] + split2[1]].ToArray();
                        }
                    }

                    // Exponential regression for starting values

                    try
                    {
                        NRe = new double[,]{{0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0},
                                                    {0}};

                        // Array resize
                        int intcount;
                        int covcount;
                        if (covariates == null)
                        {
                            covcount = 0;
                        }
                        else
                        {
                            covcount = split.Count();
                        }
                        if (int1 == null)
                        {
                            intcount = 0;
                        }
                        else
                        {
                            intcount = 1;
                        }

                        newArray = new double[covcount + intcount + 2, 1];
                        Array.Copy(NRe, newArray, newArray.Length);
                        NRe = newArray;

                        da = new double[dos.Rows.Count];
                        db = new double[dos.Rows.Count];
                        dc = new double[dos.Rows.Count];
                        dd = new double[dos.Rows.Count];
                        de = new double[dos.Rows.Count];
                        df = new double[dos.Rows.Count];
                        dg = new double[dos.Rows.Count];
                        dh = new double[dos.Rows.Count];
                        di = new double[dos.Rows.Count];
                        dj = new double[dos.Rows.Count];
                        dk = new double[dos.Rows.Count];
                        dl = new double[dos.Rows.Count];

                        d2a = new double[dos.Rows.Count];
                        d2b = new double[dos.Rows.Count];
                        d2c = new double[dos.Rows.Count];
                        d2d = new double[dos.Rows.Count];
                        d2e = new double[dos.Rows.Count];
                        d2f = new double[dos.Rows.Count];
                        d2g = new double[dos.Rows.Count];
                        d2h = new double[dos.Rows.Count];
                        d2i = new double[dos.Rows.Count];
                        d2j = new double[dos.Rows.Count];
                        d2k = new double[dos.Rows.Count];
                        d2l = new double[dos.Rows.Count];

                        dab = new double[dos.Rows.Count];
                        dac = new double[dos.Rows.Count];
                        dad = new double[dos.Rows.Count];
                        dae = new double[dos.Rows.Count];
                        daf = new double[dos.Rows.Count];
                        dag = new double[dos.Rows.Count];
                        dah = new double[dos.Rows.Count];
                        dai = new double[dos.Rows.Count];
                        daj = new double[dos.Rows.Count];
                        dak = new double[dos.Rows.Count];
                        dal = new double[dos.Rows.Count];

                        dbc = new double[dos.Rows.Count];
                        dbd = new double[dos.Rows.Count];
                        dbe = new double[dos.Rows.Count];
                        dbf = new double[dos.Rows.Count];
                        dbg = new double[dos.Rows.Count];
                        dbh = new double[dos.Rows.Count];
                        dbi = new double[dos.Rows.Count];
                        dbj = new double[dos.Rows.Count];
                        dbk = new double[dos.Rows.Count];
                        dbl = new double[dos.Rows.Count];

                        dcd = new double[dos.Rows.Count];
                        dce = new double[dos.Rows.Count];
                        dcf = new double[dos.Rows.Count];
                        dcg = new double[dos.Rows.Count];
                        dch = new double[dos.Rows.Count];
                        dci = new double[dos.Rows.Count];
                        dcj = new double[dos.Rows.Count];
                        dck = new double[dos.Rows.Count];
                        dcl = new double[dos.Rows.Count];

                        dde = new double[dos.Rows.Count];
                        ddf = new double[dos.Rows.Count];
                        ddg = new double[dos.Rows.Count];
                        ddh = new double[dos.Rows.Count];
                        ddi = new double[dos.Rows.Count];
                        ddj = new double[dos.Rows.Count];
                        ddk = new double[dos.Rows.Count];
                        ddl = new double[dos.Rows.Count];

                        def = new double[dos.Rows.Count];
                        deg = new double[dos.Rows.Count];
                        deh = new double[dos.Rows.Count];
                        dei = new double[dos.Rows.Count];
                        dej = new double[dos.Rows.Count];
                        dek = new double[dos.Rows.Count];
                        del = new double[dos.Rows.Count];

                        dfg = new double[dos.Rows.Count];
                        dfh = new double[dos.Rows.Count];
                        dfi = new double[dos.Rows.Count];
                        dfj = new double[dos.Rows.Count];
                        dfk = new double[dos.Rows.Count];
                        dfl = new double[dos.Rows.Count];

                        dgh = new double[dos.Rows.Count];
                        dgi = new double[dos.Rows.Count];
                        dgj = new double[dos.Rows.Count];
                        dgk = new double[dos.Rows.Count];
                        dgl = new double[dos.Rows.Count];

                        dhi = new double[dos.Rows.Count];
                        dhj = new double[dos.Rows.Count];
                        dhk = new double[dos.Rows.Count];
                        dhl = new double[dos.Rows.Count];

                        dij = new double[dos.Rows.Count];
                        dik = new double[dos.Rows.Count];
                        dil = new double[dos.Rows.Count];

                        djk = new double[dos.Rows.Count];
                        djl = new double[dos.Rows.Count];

                        dkl = new double[dos.Rows.Count];

                        for (int ii = 0; ii < 100; ++ii)
                        {

                            if (ii == 0)
                            {
                                b0 = 0;
                                b1 = 0;
                                b2 = 0;
                                b3 = 0;
                                b4 = 0;
                                b5 = 0;
                                b6 = 0;
                                b7 = 0;
                                b8 = 0;
                                b9 = 0;
                                b10 = 0;
                                b11 = 0;

                                ini1 = new double[,]{{b0},
                                                     {b1},
                                                     {b2},
                                                     {b3},
                                                     {b4},
                                                     {b5},
                                                     {b6},
                                                     {b7},
                                                     {b8},
                                                     {b9},
                                                     {b10},
                                                     {b11}};

                                newArray1 = new double[covcount + intcount + 2, 1];
                                Array.Copy(ini1, newArray1, newArray1.Length);
                                ini1 = newArray1;
                            }
                            else
                            {
                                b0 = NRe[0, 0];
                                b1 = NRe[1, 0];

                                if (NRe.Length >= 3)
                                {
                                    b2 = NRe[2, 0];
                                }

                                if (NRe.Length >= 4)
                                {
                                    b3 = NRe[3, 0];
                                }
                                if (NRe.Length >= 5)
                                {
                                    b4 = NRe[4, 0];
                                }
                                if (NRe.Length >= 6)
                                {
                                    b5 = NRe[5, 0];
                                }
                                if (NRe.Length >= 7)
                                {
                                    b6 = NRe[6, 0];
                                }
                                if (NRe.Length >= 8)
                                {
                                    b7 = NRe[7, 0];
                                }
                                if (NRe.Length >= 9)
                                {
                                    b8 = NRe[8, 0];
                                }
                                if (NRe.Length >= 10)
                                {
                                    b9 = NRe[9, 0];
                                }
                                if (NRe.Length >= 11)
                                {
                                    b10 = NRe[10, 0];
                                }
                                if (NRe.Length >= 12)
                                {
                                    b11 = NRe[11, 0];
                                }

                                ini1 = NRe;
                            }

                            for (int i = 0; i < dos.Rows.Count; ++i)
                            {
                                da[i] = (-(censor[i]) + (time[i] / (Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                db[i] = ((-(censor[i] * SNP[i])) + ((time[i] * SNP[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                d2a[i] = (-time[i] / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                d2b[i] = (-(time[i] * Math.Pow(SNP[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                dab[i] = (-(time[i] * SNP[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                if (covcount + intcount + 2 >= 3)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i])) + ((time[i] * Cov1[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2c[i] = (-(time[i] * Math.Pow(Cov1[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dac[i] = (-(time[i] * Cov1[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbc[i] = (-(time[i] * SNP[i] * Cov1[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                }
                                if (covcount + intcount + 2 >= 4)
                                {
                                    dd[i] = ((-(censor[i] * Cov2[i])) + ((time[i] * Cov2[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2d[i] = (-(time[i] * Math.Pow(Cov2[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dad[i] = (-(time[i] * Cov2[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbd[i] = (-(time[i] * SNP[i] * Cov2[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dcd[i] = (-(time[i] * Cov1[i] * Cov2[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                }
                                if (covcount + intcount + 2 >= 5)
                                {
                                    de[i] = ((-(censor[i] * Cov3[i])) + ((time[i] * Cov3[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2e[i] = (-(time[i] * Math.Pow(Cov3[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dae[i] = (-(time[i] * Cov3[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbe[i] = (-(time[i] * SNP[i] * Cov3[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dce[i] = (-(time[i] * Cov1[i] * Cov3[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dde[i] = (-(time[i] * Cov2[i] * Cov3[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                }
                                if (covcount + intcount + 2 >= 6)
                                {
                                    df[i] = ((-(censor[i] * Cov4[i])) + ((time[i] * Cov4[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2f[i] = (-(time[i] * Math.Pow(Cov4[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    daf[i] = (-(time[i] * Cov4[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbf[i] = (-(time[i] * SNP[i] * Cov4[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dcf[i] = (-(time[i] * Cov1[i] * Cov4[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    ddf[i] = (-(time[i] * Cov2[i] * Cov4[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    def[i] = (-(time[i] * Cov3[i] * Cov4[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                }
                                if (covcount + intcount + 2 >= 7)
                                {
                                    dg[i] = ((-(censor[i] * Cov5[i])) + ((time[i] * Cov5[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2g[i] = (-(time[i] * Math.Pow(Cov5[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dag[i] = (-(time[i] * Cov5[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbg[i] = (-(time[i] * SNP[i] * Cov5[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dcg[i] = (-(time[i] * Cov1[i] * Cov5[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    ddg[i] = (-(time[i] * Cov2[i] * Cov5[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    deg[i] = (-(time[i] * Cov3[i] * Cov5[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dfg[i] = (-(time[i] * Cov4[i] * Cov5[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                }
                                if (covcount + intcount + 2 >= 8)
                                {
                                    dh[i] = ((-(censor[i] * Cov6[i])) + ((time[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2h[i] = (-(time[i] * Math.Pow(Cov6[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dah[i] = (-(time[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbh[i] = (-(time[i] * SNP[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dch[i] = (-(time[i] * Cov1[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    ddh[i] = (-(time[i] * Cov2[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    deh[i] = (-(time[i] * Cov3[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dfh[i] = (-(time[i] * Cov4[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dgh[i] = (-(time[i] * Cov5[i] * Cov6[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                }
                                if (covcount + intcount + 2 >= 9)
                                {
                                    di[i] = ((-(censor[i] * Cov7[i])) + ((time[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2i[i] = (-(time[i] * Math.Pow(Cov7[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dai[i] = (-(time[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbi[i] = (-(time[i] * SNP[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dci[i] = (-(time[i] * Cov1[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    ddi[i] = (-(time[i] * Cov2[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dei[i] = (-(time[i] * Cov3[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dfi[i] = (-(time[i] * Cov4[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dgi[i] = (-(time[i] * Cov5[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dhi[i] = (-(time[i] * Cov6[i] * Cov7[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                }
                                if (covcount + intcount + 2 >= 10)
                                {
                                    dj[i] = ((-(censor[i] * Cov8[i])) + ((time[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2j[i] = (-(time[i] * Math.Pow(Cov8[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    daj[i] = (-(time[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbj[i] = (-(time[i] * SNP[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dcj[i] = (-(time[i] * Cov1[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    ddj[i] = (-(time[i] * Cov2[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dej[i] = (-(time[i] * Cov3[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dfj[i] = (-(time[i] * Cov4[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dgj[i] = (-(time[i] * Cov5[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dhj[i] = (-(time[i] * Cov6[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dij[i] = (-(time[i] * Cov7[i] * Cov8[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                }
                                if (covcount + intcount + 2 >= 11)
                                {
                                    dk[i] = ((-(censor[i] * Cov9[i])) + ((time[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                                    d2k[i] = (-(time[i] * Math.Pow(Cov9[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dak[i] = (-(time[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbk[i] = (-(time[i] * SNP[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dck[i] = (-(time[i] * Cov1[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    ddk[i] = (-(time[i] * Cov2[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dek[i] = (-(time[i] * Cov3[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dfk[i] = (-(time[i] * Cov4[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dgk[i] = (-(time[i] * Cov5[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dhk[i] = (-(time[i] * Cov6[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dik[i] = (-(time[i] * Cov7[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    djk[i] = (-(time[i] * Cov8[i] * Cov9[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));

                                }
                                if (covcount + intcount + 2 >= 12)
                                {
                                    d2l[i] = (-(time[i] * Math.Pow(Cov10[i], 2)) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dal[i] = (-(time[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dbl[i] = (-(time[i] * SNP[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dcl[i] = (-(time[i] * Cov1[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    ddl[i] = (-(time[i] * Cov2[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    del[i] = (-(time[i] * Cov3[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dfl[i] = (-(time[i] * Cov4[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dgl[i] = (-(time[i] * Cov5[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dhl[i] = (-(time[i] * Cov6[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dil[i] = (-(time[i] * Cov7[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    djl[i] = (-(time[i] * Cov8[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                    dkl[i] = (-(time[i] * Cov9[i] * Cov10[i]) / ((Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))));
                                }
                            }
                            sumda = da.Sum();
                            sumdb = db.Sum();
                            sumdc = dc.Sum();
                            sumdd = dd.Sum();
                            sumde = de.Sum();
                            sumdf = df.Sum();
                            sumdg = dg.Sum();
                            sumdh = dh.Sum();
                            sumdi = di.Sum();
                            sumdj = dj.Sum();
                            sumdk = dk.Sum();
                            sumdl = dl.Sum();

                            sumd2a = d2a.Sum();
                            sumd2b = d2b.Sum();
                            sumd2c = d2c.Sum();
                            sumd2d = d2d.Sum();
                            sumd2e = d2e.Sum();
                            sumd2f = d2f.Sum();
                            sumd2g = d2g.Sum();
                            sumd2h = d2h.Sum();
                            sumd2i = d2i.Sum();
                            sumd2j = d2j.Sum();
                            sumd2k = d2k.Sum();
                            sumd2l = d2l.Sum();

                            sumdab = dab.Sum();
                            sumdac = dac.Sum();
                            sumdad = dad.Sum();
                            sumdae = dae.Sum();
                            sumdaf = daf.Sum();
                            sumdag = dag.Sum();
                            sumdah = dah.Sum();
                            sumdai = dai.Sum();
                            sumdaj = daj.Sum();
                            sumdak = dak.Sum();
                            sumdal = dal.Sum();

                            sumdbc = dbc.Sum();
                            sumdbd = dbd.Sum();
                            sumdbe = dbe.Sum();
                            sumdbf = dbf.Sum();
                            sumdbg = dbg.Sum();
                            sumdbh = dbh.Sum();
                            sumdbi = dbi.Sum();
                            sumdbj = dbj.Sum();
                            sumdbk = dbk.Sum();
                            sumdbl = dbl.Sum();

                            sumdcd = dcd.Sum();
                            sumdce = dce.Sum();
                            sumdcf = dcf.Sum();
                            sumdcg = dcg.Sum();
                            sumdch = dch.Sum();
                            sumdci = dci.Sum();
                            sumdcj = dcj.Sum();
                            sumdck = dck.Sum();
                            sumdcl = dcl.Sum();

                            sumdde = dde.Sum();
                            sumddf = ddf.Sum();
                            sumddg = ddg.Sum();
                            sumddh = ddh.Sum();
                            sumddi = ddi.Sum();
                            sumddj = ddj.Sum();
                            sumddk = ddk.Sum();
                            sumddl = ddl.Sum();

                            sumdef = def.Sum();
                            sumdeg = deg.Sum();
                            sumdeh = deh.Sum();
                            sumdei = dei.Sum();
                            sumdej = dej.Sum();
                            sumdek = dek.Sum();
                            sumdel = del.Sum();

                            sumdfg = dfg.Sum();
                            sumdfh = dfh.Sum();
                            sumdfi = dfi.Sum();
                            sumdfj = dfj.Sum();
                            sumdfk = dfk.Sum();
                            sumdfl = dfl.Sum();

                            sumdgh = dgh.Sum();
                            sumdgi = dgi.Sum();
                            sumdgj = dgj.Sum();
                            sumdgk = dgk.Sum();
                            sumdgl = dgl.Sum();

                            sumdhi = dhi.Sum();
                            sumdhj = dhj.Sum();
                            sumdhk = dhk.Sum();
                            sumdhl = dhl.Sum();

                            sumdij = dij.Sum();
                            sumdik = dik.Sum();
                            sumdil = dil.Sum();

                            sumdjk = djk.Sum();
                            sumdjl = djl.Sum();

                            sumdkl = dkl.Sum();


                            H = new double[12, 12] {{sumd2a,sumdab,sumdac,sumdad,sumdae,sumdaf,sumdag,sumdah,sumdai,sumdaj,sumdak,sumdal},
                                                               {sumdab,sumd2b,sumdbc,sumdbd,sumdbe,sumdbf,sumdbg,sumdbh,sumdbi,sumdbj,sumdbk,sumdbl},
                                                               {sumdac,sumdbc,sumd2c,sumdcd,sumdce,sumdcf,sumdcg,sumdch,sumdci,sumdcj,sumdck,sumdcl},
                                                               {sumdad,sumdbd,sumdcd,sumd2d,sumdde,sumddf,sumddg,sumddh,sumddi,sumddj,sumddk,sumddl},
                                                               {sumdae,sumdbe,sumdce,sumdde,sumd2e,sumdef,sumdeg,sumdeh,sumdei,sumdej,sumdek,sumdel},
                                                               {sumdaf,sumdbf,sumdcf,sumddf,sumdef,sumd2f,sumdfg,sumdfh,sumdfi,sumdfj,sumdfk,sumdfl},
                                                               {sumdag,sumdbg,sumdcg,sumddg,sumdeg,sumdfg,sumd2g,sumdgh,sumdgi,sumdgj,sumdgk,sumdgl},
                                                               {sumdah,sumdbh,sumdch,sumddh,sumdeh,sumdfh,sumdgh,sumd2h,sumdhi,sumdhj,sumdhk,sumdhl},
                                                               {sumdai,sumdbi,sumdci,sumddi,sumdei,sumdfi,sumdgi,sumdhi,sumd2i,sumdij,sumdik,sumdil},
                                                               {sumdaj,sumdbj,sumdcj,sumddj,sumdej,sumdfj,sumdgj,sumdhj,sumdij,sumd2j,sumdjk,sumdjl},
                                                               {sumdak,sumdbk,sumdck,sumddk,sumdek,sumdfk,sumdgk,sumdhk,sumdik,sumdjk,sumd2k,sumdkl},
                                                               {sumdal,sumdbl,sumdcl,sumddl,sumdel,sumdfl,sumdgl,sumdhl,sumdil,sumdjl,sumdkl,sumd2l}};

                            newArray2 = new double[covcount + intcount + 2, covcount + intcount + 2];
                            newArray2[0, 0] = H[0, 0];
                            newArray2[0, 1] = H[0, 1];
                            newArray2[1, 0] = H[1, 0];
                            newArray2[1, 1] = H[1, 1];

                            if (covcount + intcount + 2 >= 3)
                            {
                                newArray2[0, 2] = H[0, 2];
                                newArray2[1, 2] = H[1, 2];
                                newArray2[2, 2] = H[2, 2];
                                newArray2[2, 0] = H[2, 0];
                                newArray2[2, 1] = H[2, 1];
                            }
                            if (covcount + intcount + 2 >= 4)
                            {
                                newArray2[0, 3] = H[0, 3];
                                newArray2[1, 3] = H[1, 3];
                                newArray2[2, 3] = H[2, 3];
                                newArray2[3, 3] = H[3, 3];
                                newArray2[3, 0] = H[3, 0];
                                newArray2[3, 1] = H[3, 1];
                                newArray2[3, 2] = H[3, 2];
                            }
                            if (covcount + intcount + 2 >= 5)
                            {
                                newArray2[0, 4] = H[0, 4];
                                newArray2[1, 4] = H[1, 4];
                                newArray2[2, 4] = H[2, 4];
                                newArray2[3, 4] = H[3, 4];
                                newArray2[4, 4] = H[4, 4];
                                newArray2[4, 0] = H[4, 0];
                                newArray2[4, 1] = H[4, 1];
                                newArray2[4, 2] = H[4, 2];
                                newArray2[4, 3] = H[4, 3];
                            }
                            if (covcount + intcount + 2 >= 6)
                            {
                                newArray2[0, 5] = H[0, 5];
                                newArray2[1, 5] = H[1, 5];
                                newArray2[2, 5] = H[2, 5];
                                newArray2[3, 5] = H[3, 5];
                                newArray2[4, 5] = H[4, 5];
                                newArray2[5, 5] = H[5, 5];
                                newArray2[5, 0] = H[5, 0];
                                newArray2[5, 1] = H[5, 1];
                                newArray2[5, 2] = H[5, 2];
                                newArray2[5, 3] = H[5, 3];
                                newArray2[5, 4] = H[5, 4];
                            }
                            if (covcount + intcount + 2 >= 7)
                            {
                                newArray2[0, 6] = H[0, 6];
                                newArray2[1, 6] = H[1, 6];
                                newArray2[2, 6] = H[2, 6];
                                newArray2[3, 6] = H[3, 6];
                                newArray2[4, 6] = H[4, 6];
                                newArray2[5, 6] = H[5, 6];
                                newArray2[6, 6] = H[6, 6];
                                newArray2[6, 0] = H[6, 0];
                                newArray2[6, 1] = H[6, 1];
                                newArray2[6, 2] = H[6, 2];
                                newArray2[6, 3] = H[6, 3];
                                newArray2[6, 4] = H[6, 4];
                                newArray2[6, 5] = H[6, 5];
                            }
                            if (covcount + intcount + 2 >= 8)
                            {
                                newArray2[0, 7] = H[0, 7];
                                newArray2[1, 7] = H[1, 7];
                                newArray2[2, 7] = H[2, 7];
                                newArray2[3, 7] = H[3, 7];
                                newArray2[4, 7] = H[4, 7];
                                newArray2[5, 7] = H[5, 7];
                                newArray2[6, 7] = H[6, 7];
                                newArray2[7, 7] = H[7, 7];
                                newArray2[7, 0] = H[7, 0];
                                newArray2[7, 1] = H[7, 1];
                                newArray2[7, 2] = H[7, 2];
                                newArray2[7, 3] = H[7, 3];
                                newArray2[7, 4] = H[7, 4];
                                newArray2[7, 5] = H[7, 5];
                                newArray2[7, 6] = H[7, 6];
                            }
                            if (covcount + intcount + 2 >= 9)
                            {
                                newArray2[0, 8] = H[0, 8];
                                newArray2[1, 8] = H[1, 8];
                                newArray2[2, 8] = H[2, 8];
                                newArray2[3, 8] = H[3, 8];
                                newArray2[4, 8] = H[4, 8];
                                newArray2[5, 8] = H[5, 8];
                                newArray2[6, 8] = H[6, 8];
                                newArray2[7, 8] = H[7, 8];
                                newArray2[8, 8] = H[8, 8];
                                newArray2[8, 0] = H[8, 0];
                                newArray2[8, 1] = H[8, 1];
                                newArray2[8, 2] = H[8, 2];
                                newArray2[8, 3] = H[8, 3];
                                newArray2[8, 4] = H[8, 4];
                                newArray2[8, 5] = H[8, 5];
                                newArray2[8, 6] = H[8, 6];
                                newArray2[8, 7] = H[8, 7];
                            }
                            if (covcount + intcount + 2 >= 10)
                            {
                                newArray2[0, 9] = H[0, 9];
                                newArray2[1, 9] = H[1, 9];
                                newArray2[2, 9] = H[2, 9];
                                newArray2[3, 9] = H[3, 9];
                                newArray2[4, 9] = H[4, 9];
                                newArray2[5, 9] = H[5, 9];
                                newArray2[6, 9] = H[6, 9];
                                newArray2[7, 9] = H[7, 9];
                                newArray2[8, 9] = H[8, 9];
                                newArray2[9, 9] = H[9, 9];
                                newArray2[9, 0] = H[9, 0];
                                newArray2[9, 1] = H[9, 1];
                                newArray2[9, 2] = H[9, 2];
                                newArray2[9, 3] = H[9, 3];
                                newArray2[9, 4] = H[9, 4];
                                newArray2[9, 5] = H[9, 5];
                                newArray2[9, 6] = H[9, 6];
                                newArray2[9, 7] = H[9, 7];
                                newArray2[9, 8] = H[9, 8];
                            }
                            if (covcount + intcount + 2 >= 11)
                            {
                                newArray2[0, 10] = H[0, 10];
                                newArray2[1, 10] = H[1, 10];
                                newArray2[2, 10] = H[2, 10];
                                newArray2[3, 10] = H[3, 10];
                                newArray2[4, 10] = H[4, 10];
                                newArray2[5, 10] = H[5, 10];
                                newArray2[6, 10] = H[6, 10];
                                newArray2[7, 10] = H[7, 10];
                                newArray2[8, 10] = H[8, 10];
                                newArray2[9, 10] = H[9, 10];
                                newArray2[10, 10] = H[10, 10];
                                newArray2[10, 0] = H[10, 0];
                                newArray2[10, 1] = H[10, 1];
                                newArray2[10, 2] = H[10, 2];
                                newArray2[10, 3] = H[10, 3];
                                newArray2[10, 4] = H[10, 4];
                                newArray2[10, 5] = H[10, 5];
                                newArray2[10, 6] = H[10, 6];
                                newArray2[10, 7] = H[10, 7];
                                newArray2[10, 8] = H[10, 8];
                                newArray2[10, 9] = H[10, 9];
                            }
                            if (covcount + intcount + 2 >= 12)
                            {
                                newArray2[0, 11] = H[0, 11];
                                newArray2[1, 11] = H[1, 11];
                                newArray2[2, 11] = H[2, 11];
                                newArray2[3, 11] = H[3, 11];
                                newArray2[4, 11] = H[4, 11];
                                newArray2[5, 11] = H[5, 11];
                                newArray2[6, 11] = H[6, 11];
                                newArray2[7, 11] = H[7, 11];
                                newArray2[8, 11] = H[8, 11];
                                newArray2[9, 11] = H[9, 11];
                                newArray2[10, 11] = H[10, 11];
                                newArray2[11, 0] = H[11, 0];
                                newArray2[11, 1] = H[11, 1];
                                newArray2[11, 2] = H[11, 2];
                                newArray2[11, 3] = H[11, 3];
                                newArray2[11, 4] = H[11, 4];
                                newArray2[11, 5] = H[11, 5];
                                newArray2[11, 6] = H[11, 6];
                                newArray2[11, 7] = H[11, 7];
                                newArray2[11, 8] = H[11, 8];
                                newArray2[11, 9] = H[11, 9];
                                newArray2[11, 10] = H[11, 10];
                                newArray2[11, 11] = H[11, 11];
                            }

                            H = newArray2;

                            IH = H.Inverse();

                            Z = new double[,]{{sumda},
                                                       {sumdb},
                                                       {sumdc},
                                                       {sumdd},
                                                       {sumde},
                                                       {sumdf},
                                                       {sumdg},
                                                       {sumdh},
                                                       {sumdi},
                                                       {sumdj},
                                                       {sumdk},
                                                       {sumdl}};

                            newArray3 = new double[covcount + intcount + 2, 1];
                            Array.Copy(Z, newArray3, newArray3.Length);
                            Z = newArray3;

                            // newton raphson

                            NRe = ini1.Subtract(IH.Multiply(Z));
                            if (NRe[0, 0].Equals(Double.NaN) || NRe[1, 0].Equals(Double.NaN))
                            {
                                break;
                            }
                            if (Math.Abs(NRe[1, 0] - b1) < 0.000001)
                            { break; }
                            else
                            {
                                continue;
                            }

                        }

                        if (NRe[0, 0].Equals(Double.NaN) || NRe[1, 0].Equals(Double.NaN))
                        {
                            Console.WriteLine("NaN's produced, " + dos.Columns[0] + " analysis did not converge");
                            goto SKIP2;
                        }
                        NR = new double[covcount + intcount + 3, 1];
                        ini = new double[covcount + intcount + 3, 1];

                        if (covcount + intcount + 3 <= 3)
                        {

                            ini = new double[,]{{NRe[0, 0]},
                                             {NRe[1, 0]},
                                             {0.1}};
                        }
                        if (covcount + intcount + 3 == 4)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 5)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 6)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 7)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 8)
                        {
                            ini = new double[,] {{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 9)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 10)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 11)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {NRe[9, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 12)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {NRe[9, 0]},
                                                    {NRe[10, 0]},
                                                    {0.1}};
                        }
                        if (covcount + intcount + 3 == 13)
                        {
                            ini = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {NRe[9, 0]},
                                                    {NRe[10, 0]},
                                                    {NRe[11, 0]},
                                                    {0.1}};
                        }

                        dm = new double[dos.Rows.Count];
                        d2m = new double[dos.Rows.Count];
                        dma = new double[dos.Rows.Count];
                        dmb = new double[dos.Rows.Count];
                        dmc = new double[dos.Rows.Count];
                        dmd = new double[dos.Rows.Count];
                        dme = new double[dos.Rows.Count];
                        dmf = new double[dos.Rows.Count];
                        dmg = new double[dos.Rows.Count];
                        dmh = new double[dos.Rows.Count];
                        dmi = new double[dos.Rows.Count];
                        dmj = new double[dos.Rows.Count];
                        dmk = new double[dos.Rows.Count];
                        dml = new double[dos.Rows.Count];

                        double sig1 = 0.1;

                        for (int ii = 0; ii < 100000; ++ii)
                        {

                            if (ii == 0)
                            {
                                b0 = NRe[0, 0];
                                b1 = NRe[1, 0];
                                sig = 0.1;

                                try
                                {
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                    b6 = NRe[6, 0];
                                    b7 = NRe[7, 0];
                                    b8 = NRe[8, 0];
                                    b9 = NRe[9, 0];
                                    b10 = NRe[10, 0];
                                    b11 = NRe[11, 0];

                                }
                                catch (Exception betas)
                                {
                                }
                            }
                            else
                            {
                                ini = NR;

                                b0 = NR[0, 0];
                                b1 = NR[1, 0];

                                if (NR.Length <= 3)
                                {
                                    sig = NR[2, 0];
                                }

                                if (NR.Length == 4)
                                {
                                    b2 = NR[2, 0];
                                    sig = NR[3, 0];
                                }

                                if (NR.Length == 5)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[4, 0];
                                }
                                if (NR.Length == 6)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[5, 0];
                                    b4 = NR[4, 0];
                                }
                                if (NR.Length == 7)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[6, 0];
                                    b4 = NR[4, 0];
                                    b5 = NR[5, 0];
                                }
                                if (NR.Length == 8)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[7, 0];
                                    b4 = NR[4, 0];
                                    b5 = NR[5, 0];
                                    b6 = NR[6, 0];
                                }
                                if (NR.Length == 9)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[8, 0];
                                    b4 = NR[4, 0];
                                    b5 = NR[5, 0];
                                    b6 = NR[6, 0];
                                    b7 = NR[7, 0];
                                }
                                if (NR.Length == 10)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[9, 0];
                                    b4 = NR[4, 0];
                                    b5 = NR[5, 0];
                                    b6 = NR[6, 0];
                                    b7 = NR[7, 0];
                                    b8 = NR[8, 0];
                                }
                                if (NR.Length == 11)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[10, 0];
                                    b4 = NR[4, 0];
                                    b5 = NR[5, 0];
                                    b6 = NR[6, 0];
                                    b7 = NR[7, 0];
                                    b8 = NR[8, 0];
                                    b9 = NR[9, 0];
                                }
                                if (NR.Length == 12)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[11, 0];
                                    b4 = NR[4, 0];
                                    b5 = NR[5, 0];
                                    b6 = NR[6, 0];
                                    b7 = NR[7, 0];
                                    b8 = NR[8, 0];
                                    b9 = NR[9, 0];
                                    b10 = NR[10, 0];
                                }
                                if (NR.Length == 13)
                                {
                                    b2 = NR[2, 0];
                                    b3 = NR[3, 0];
                                    sig = NR[12, 0];
                                    b4 = NR[4, 0];
                                    b5 = NR[5, 0];
                                    b6 = NR[6, 0];
                                    b7 = NR[7, 0];
                                    b8 = NR[8, 0];
                                    b9 = NR[9, 0];
                                    b10 = NR[10, 0];
                                    b11 = NR[11, 0];
                                }
                            }

                        Try:;
                            if (valid == 1)
                            {
                                sig = 16;
                            }
                            if ((sig >= 15 || sig < 0) || (NR[0, 0] >= 30 || NR[0, 0] <= -30) || (NR[1, 0] >= 30 || NR[1, 0] <= -30))

                            {
                                itercount++;
                                if (itercount >= 80)
                                {
                                    Console.WriteLine("NaN's produced, " + dos.Columns[0] + " analysis did not converge");
                                    goto SKIP2;
                                }

                                sig1 += 0.1;
                                sig = sig1;

                                if (covcount + intcount + 3 <= 3)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {sig}};

                                    ini = ini2;
                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                }
                                if (covcount + intcount + 3 == 4)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {sig}};

                                    ini = ini2;
                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                }
                                if (covcount + intcount + 3 == 5)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {sig}};

                                    ini = ini2;
                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                }
                                if (covcount + intcount + 3 == 6)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {sig}};

                                    ini = ini2;
                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                }
                                if (covcount + intcount + 3 == 7)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {sig}};

                                    ini = ini2;
                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                }
                                if (covcount + intcount + 3 == 8)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {sig}};

                                    ini = ini2;
                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                    b6 = NRe[6, 0];
                                }
                                if (covcount + intcount + 3 == 9)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {sig}};

                                    ini = ini2;
                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                    b6 = NRe[6, 0];
                                    b7 = NRe[7, 0];
                                }
                                if (covcount + intcount + 3 == 10)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {sig}};

                                    ini = ini2;

                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                    b6 = NRe[6, 0];
                                    b7 = NRe[7, 0];
                                    b8 = NRe[8, 0];
                                }
                                if (covcount + intcount + 3 == 11)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {NRe[9, 0]},
                                                    {sig}};

                                    ini = ini2;

                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                    b6 = NRe[6, 0];
                                    b7 = NRe[7, 0];
                                    b8 = NRe[8, 0];
                                    b9 = NRe[9, 0];
                                }
                                if (covcount + intcount + 3 == 12)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {NRe[9, 0]},
                                                    {NRe[10, 0]},
                                                    {sig}};

                                    ini = ini2;

                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                    b6 = NRe[6, 0];
                                    b7 = NRe[7, 0];
                                    b8 = NRe[8, 0];
                                    b9 = NRe[9, 0];
                                    b10 = NRe[10, 0];
                                }
                                if (covcount + intcount + 3 == 13)
                                {
                                    ini2 = new double[,]{{NRe[0, 0]},
                                                    {NRe[1, 0]},
                                                    {NRe[2, 0]},
                                                    {NRe[3, 0]},
                                                    {NRe[4, 0]},
                                                    {NRe[5, 0]},
                                                    {NRe[6, 0]},
                                                    {NRe[7, 0]},
                                                    {NRe[8, 0]},
                                                    {NRe[9, 0]},
                                                    {NRe[10, 0]},
                                                    {NRe[11, 0]},
                                                    {sig}};

                                    ini = ini2;

                                    b0 = NRe[0, 0];
                                    b1 = NRe[1, 0];
                                    b2 = NRe[2, 0];
                                    b3 = NRe[3, 0];
                                    b4 = NRe[4, 0];
                                    b5 = NRe[5, 0];
                                    b6 = NRe[6, 0];
                                    b7 = NRe[7, 0];
                                    b8 = NRe[8, 0];
                                    b9 = NRe[9, 0];
                                    b10 = NRe[10, 0];
                                    b11 = NRe[11, 0];
                                }
                            }



                            for (int i = 0; i < dos.Rows.Count; ++i)
                            {
                                da[i] = (-(censor[i] / sig) + ((Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                db[i] = ((-(censor[i] * SNP[i]) / sig) + ((SNP[i] * Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                d2a[i] = ((-Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                d2b[i] = ((-(Math.Pow(SNP[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)))) / (Math.Pow(sig, 2)); //check these brackets
                                dab[i] = (-SNP[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));

                                if (covcount + intcount + 3 == 3)
                                {
                                    // sigma
                                    dc[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2c[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dac[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbc[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                }

                                if (covcount + intcount + 3 == 4)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2d[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dad[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbd[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dcd[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                }

                                if (covcount + intcount + 3 == 5)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2e[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dae[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbe[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dce[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dde[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }

                                if (covcount + intcount + 3 == 6)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2f[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    daf[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbf[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dcf[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    ddf[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    def[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }

                                if (covcount + intcount + 3 == 7)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-(censor[i] * Cov4[i]) / sig) + (Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2f[i] = ((-(Math.Pow(Cov4[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daf[i] = (-Cov4[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbf[i] = (-SNP[i] * Cov4[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcf[i] = (-Cov1[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddf[i] = (-Cov2[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    def[i] = (-Cov3[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dg[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2g[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dag[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbg[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dcg[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    ddg[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    deg[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dfg[i] = ((censor[i] * Cov4[i]) / (Math.Pow(sig, 2))) - ((Cov4[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }

                                if (covcount + intcount + 3 == 8)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-(censor[i] * Cov4[i]) / sig) + (Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2f[i] = ((-(Math.Pow(Cov4[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daf[i] = (-Cov4[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbf[i] = (-SNP[i] * Cov4[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcf[i] = (-Cov1[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddf[i] = (-Cov2[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    def[i] = (-Cov3[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dg[i] = ((-(censor[i] * Cov5[i]) / sig) + (Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2g[i] = ((-(Math.Pow(Cov5[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dag[i] = (-Cov5[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbg[i] = (-SNP[i] * Cov5[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcg[i] = (-Cov1[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddg[i] = (-Cov2[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deg[i] = (-Cov3[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfg[i] = (-Cov4[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dh[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2h[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dah[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbh[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dch[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    ddh[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    deh[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dfh[i] = ((censor[i] * Cov4[i]) / (Math.Pow(sig, 2))) - ((Cov4[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dgh[i] = ((censor[i] * Cov5[i]) / (Math.Pow(sig, 2))) - ((Cov5[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }

                                if (covcount + intcount + 3 == 9)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-(censor[i] * Cov4[i]) / sig) + (Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2f[i] = ((-(Math.Pow(Cov4[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daf[i] = (-Cov4[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbf[i] = (-SNP[i] * Cov4[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcf[i] = (-Cov1[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddf[i] = (-Cov2[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    def[i] = (-Cov3[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dg[i] = ((-(censor[i] * Cov5[i]) / sig) + (Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2g[i] = ((-(Math.Pow(Cov5[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dag[i] = (-Cov5[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbg[i] = (-SNP[i] * Cov5[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcg[i] = (-Cov1[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddg[i] = (-Cov2[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deg[i] = (-Cov3[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfg[i] = (-Cov4[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dh[i] = ((-(censor[i] * Cov6[i]) / sig) + (Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2h[i] = ((-(Math.Pow(Cov6[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dah[i] = (-Cov6[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbh[i] = (-SNP[i] * Cov6[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dch[i] = (-Cov1[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddh[i] = (-Cov2[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deh[i] = (-Cov3[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfh[i] = (-Cov4[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgh[i] = (-Cov5[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    di[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2i[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dai[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbi[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dci[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    ddi[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dei[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dfi[i] = ((censor[i] * Cov4[i]) / (Math.Pow(sig, 2))) - ((Cov4[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dgi[i] = ((censor[i] * Cov5[i]) / (Math.Pow(sig, 2))) - ((Cov5[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dhi[i] = ((censor[i] * Cov6[i]) / (Math.Pow(sig, 2))) - ((Cov6[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }
                                if (covcount + intcount + 3 == 10)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-(censor[i] * Cov4[i]) / sig) + (Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2f[i] = ((-(Math.Pow(Cov4[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daf[i] = (-Cov4[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbf[i] = (-SNP[i] * Cov4[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcf[i] = (-Cov1[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddf[i] = (-Cov2[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    def[i] = (-Cov3[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dg[i] = ((-(censor[i] * Cov5[i]) / sig) + (Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2g[i] = ((-(Math.Pow(Cov5[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dag[i] = (-Cov5[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbg[i] = (-SNP[i] * Cov5[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcg[i] = (-Cov1[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddg[i] = (-Cov2[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deg[i] = (-Cov3[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfg[i] = (-Cov4[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dh[i] = ((-(censor[i] * Cov6[i]) / sig) + (Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2h[i] = ((-(Math.Pow(Cov6[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dah[i] = (-Cov6[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbh[i] = (-SNP[i] * Cov6[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dch[i] = (-Cov1[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddh[i] = (-Cov2[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deh[i] = (-Cov3[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfh[i] = (-Cov4[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgh[i] = (-Cov5[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    di[i] = ((-(censor[i] * Cov7[i]) / sig) + (Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2i[i] = ((-(Math.Pow(Cov7[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dai[i] = (-Cov7[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbi[i] = (-SNP[i] * Cov7[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dci[i] = (-Cov1[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddi[i] = (-Cov2[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dei[i] = (-Cov3[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfi[i] = (-Cov4[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgi[i] = (-Cov5[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhi[i] = (-Cov6[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dj[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2j[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    daj[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbj[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dcj[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    ddj[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dej[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dfj[i] = ((censor[i] * Cov4[i]) / (Math.Pow(sig, 2))) - ((Cov4[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dgj[i] = ((censor[i] * Cov5[i]) / (Math.Pow(sig, 2))) - ((Cov5[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dhj[i] = ((censor[i] * Cov6[i]) / (Math.Pow(sig, 2))) - ((Cov6[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dij[i] = ((censor[i] * Cov7[i]) / (Math.Pow(sig, 2))) - ((Cov7[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }
                                if (covcount + intcount + 3 == 11)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-(censor[i] * Cov4[i]) / sig) + (Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2f[i] = ((-(Math.Pow(Cov4[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daf[i] = (-Cov4[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbf[i] = (-SNP[i] * Cov4[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcf[i] = (-Cov1[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddf[i] = (-Cov2[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    def[i] = (-Cov3[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dg[i] = ((-(censor[i] * Cov5[i]) / sig) + (Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2g[i] = ((-(Math.Pow(Cov5[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dag[i] = (-Cov5[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbg[i] = (-SNP[i] * Cov5[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcg[i] = (-Cov1[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddg[i] = (-Cov2[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deg[i] = (-Cov3[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfg[i] = (-Cov4[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dh[i] = ((-(censor[i] * Cov6[i]) / sig) + (Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2h[i] = ((-(Math.Pow(Cov6[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dah[i] = (-Cov6[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbh[i] = (-SNP[i] * Cov6[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dch[i] = (-Cov1[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddh[i] = (-Cov2[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deh[i] = (-Cov3[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfh[i] = (-Cov4[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgh[i] = (-Cov5[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    di[i] = ((-(censor[i] * Cov7[i]) / sig) + (Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2i[i] = ((-(Math.Pow(Cov7[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dai[i] = (-Cov7[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbi[i] = (-SNP[i] * Cov7[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dci[i] = (-Cov1[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddi[i] = (-Cov2[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dei[i] = (-Cov3[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfi[i] = (-Cov4[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgi[i] = (-Cov5[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhi[i] = (-Cov6[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dj[i] = ((-(censor[i] * Cov8[i]) / sig) + (Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2j[i] = ((-(Math.Pow(Cov8[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daj[i] = (-Cov8[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbj[i] = (-SNP[i] * Cov8[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcj[i] = (-Cov1[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddj[i] = (-Cov2[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dej[i] = (-Cov3[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfj[i] = (-Cov4[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgj[i] = (-Cov5[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhj[i] = (-Cov6[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dij[i] = (-Cov7[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dk[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2k[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dak[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbk[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dck[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    ddk[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dek[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dfk[i] = ((censor[i] * Cov4[i]) / (Math.Pow(sig, 2))) - ((Cov4[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dgk[i] = ((censor[i] * Cov5[i]) / (Math.Pow(sig, 2))) - ((Cov5[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dhk[i] = ((censor[i] * Cov6[i]) / (Math.Pow(sig, 2))) - ((Cov6[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dik[i] = ((censor[i] * Cov7[i]) / (Math.Pow(sig, 2))) - ((Cov7[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    djk[i] = ((censor[i] * Cov8[i]) / (Math.Pow(sig, 2))) - ((Cov8[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }

                                if (covcount + intcount + 3 == 12)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-(censor[i] * Cov4[i]) / sig) + (Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2f[i] = ((-(Math.Pow(Cov4[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daf[i] = (-Cov4[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbf[i] = (-SNP[i] * Cov4[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcf[i] = (-Cov1[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddf[i] = (-Cov2[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    def[i] = (-Cov3[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dg[i] = ((-(censor[i] * Cov5[i]) / sig) + (Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2g[i] = ((-(Math.Pow(Cov5[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dag[i] = (-Cov5[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbg[i] = (-SNP[i] * Cov5[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcg[i] = (-Cov1[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddg[i] = (-Cov2[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deg[i] = (-Cov3[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfg[i] = (-Cov4[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dh[i] = ((-(censor[i] * Cov6[i]) / sig) + (Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2h[i] = ((-(Math.Pow(Cov6[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dah[i] = (-Cov6[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbh[i] = (-SNP[i] * Cov6[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dch[i] = (-Cov1[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddh[i] = (-Cov2[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deh[i] = (-Cov3[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfh[i] = (-Cov4[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgh[i] = (-Cov5[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    di[i] = ((-(censor[i] * Cov7[i]) / sig) + (Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2i[i] = ((-(Math.Pow(Cov7[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dai[i] = (-Cov7[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbi[i] = (-SNP[i] * Cov7[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dci[i] = (-Cov1[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddi[i] = (-Cov2[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dei[i] = (-Cov3[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfi[i] = (-Cov4[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgi[i] = (-Cov5[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhi[i] = (-Cov6[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dj[i] = ((-(censor[i] * Cov8[i]) / sig) + (Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2j[i] = ((-(Math.Pow(Cov8[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daj[i] = (-Cov8[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbj[i] = (-SNP[i] * Cov8[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcj[i] = (-Cov1[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddj[i] = (-Cov2[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dej[i] = (-Cov3[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfj[i] = (-Cov4[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgj[i] = (-Cov5[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhj[i] = (-Cov6[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dij[i] = (-Cov7[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dk[i] = ((-(censor[i] * Cov9[i]) / sig) + (Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2k[i] = ((-(Math.Pow(Cov9[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dak[i] = (-Cov9[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbk[i] = (-SNP[i] * Cov9[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dck[i] = (-Cov1[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddk[i] = (-Cov2[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dek[i] = (-Cov3[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfk[i] = (-Cov4[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgk[i] = (-Cov5[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhk[i] = (-Cov6[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dik[i] = (-Cov7[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    djk[i] = (-Cov8[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dl[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2l[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dal[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dbl[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dcl[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    ddl[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    del[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dfl[i] = ((censor[i] * Cov4[i]) / (Math.Pow(sig, 2))) - ((Cov4[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dgl[i] = ((censor[i] * Cov5[i]) / (Math.Pow(sig, 2))) - ((Cov5[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dhl[i] = ((censor[i] * Cov6[i]) / (Math.Pow(sig, 2))) - ((Cov6[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dil[i] = ((censor[i] * Cov7[i]) / (Math.Pow(sig, 2))) - ((Cov7[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    djl[i] = ((censor[i] * Cov8[i]) / (Math.Pow(sig, 2))) - ((Cov8[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dkl[i] = ((censor[i] * Cov9[i]) / (Math.Pow(sig, 2))) - ((Cov9[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));

                                }
                                if (covcount + intcount + 3 == 13)
                                {
                                    dc[i] = ((-(censor[i] * Cov1[i]) / sig) + (Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2c[i] = ((-(Math.Pow(Cov1[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dac[i] = (-Cov1[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbc[i] = (-SNP[i] * Cov1[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dd[i] = ((-(censor[i] * Cov2[i]) / sig) + (Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2d[i] = ((-(Math.Pow(Cov2[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dad[i] = (-Cov2[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbd[i] = (-SNP[i] * Cov2[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcd[i] = (-Cov1[i] * Cov2[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    de[i] = ((-(censor[i] * Cov3[i]) / sig) + (Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2e[i] = ((-(Math.Pow(Cov3[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dae[i] = (-Cov3[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbe[i] = (-SNP[i] * Cov3[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dce[i] = (-Cov1[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dde[i] = (-Cov2[i] * Cov3[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    df[i] = ((-(censor[i] * Cov4[i]) / sig) + (Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2f[i] = ((-(Math.Pow(Cov4[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daf[i] = (-Cov4[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbf[i] = (-SNP[i] * Cov4[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcf[i] = (-Cov1[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddf[i] = (-Cov2[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    def[i] = (-Cov3[i] * Cov4[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dg[i] = ((-(censor[i] * Cov5[i]) / sig) + (Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2g[i] = ((-(Math.Pow(Cov5[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dag[i] = (-Cov5[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbg[i] = (-SNP[i] * Cov5[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcg[i] = (-Cov1[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddg[i] = (-Cov2[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deg[i] = (-Cov3[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfg[i] = (-Cov4[i] * Cov5[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dh[i] = ((-(censor[i] * Cov6[i]) / sig) + (Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2h[i] = ((-(Math.Pow(Cov6[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dah[i] = (-Cov6[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbh[i] = (-SNP[i] * Cov6[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dch[i] = (-Cov1[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddh[i] = (-Cov2[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    deh[i] = (-Cov3[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfh[i] = (-Cov4[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgh[i] = (-Cov5[i] * Cov6[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    di[i] = ((-(censor[i] * Cov7[i]) / sig) + (Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2i[i] = ((-(Math.Pow(Cov7[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dai[i] = (-Cov7[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbi[i] = (-SNP[i] * Cov7[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dci[i] = (-Cov1[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddi[i] = (-Cov2[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dei[i] = (-Cov3[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfi[i] = (-Cov4[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgi[i] = (-Cov5[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhi[i] = (-Cov6[i] * Cov7[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dj[i] = ((-(censor[i] * Cov8[i]) / sig) + (Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2j[i] = ((-(Math.Pow(Cov8[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    daj[i] = (-Cov8[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbj[i] = (-SNP[i] * Cov8[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcj[i] = (-Cov1[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddj[i] = (-Cov2[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dej[i] = (-Cov3[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfj[i] = (-Cov4[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgj[i] = (-Cov5[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhj[i] = (-Cov6[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dij[i] = (-Cov7[i] * Cov8[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dk[i] = ((-(censor[i] * Cov9[i]) / sig) + (Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2k[i] = ((-(Math.Pow(Cov9[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dak[i] = (-Cov9[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbk[i] = (-SNP[i] * Cov9[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dck[i] = (-Cov1[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddk[i] = (-Cov2[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dek[i] = (-Cov3[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfk[i] = (-Cov4[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgk[i] = (-Cov5[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhk[i] = (-Cov6[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dik[i] = (-Cov7[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    djk[i] = (-Cov8[i] * Cov9[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dl[i] = ((-(censor[i] * Cov10[i]) / sig) + (Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / sig);
                                    d2l[i] = ((-(Math.Pow(Cov10[i], 2)) * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2)));
                                    dal[i] = (-Cov10[i] * (Math.Pow(time[i], (1 / sig))) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2));
                                    dbl[i] = (-SNP[i] * Cov10[i] * (Math.Pow(time[i], (1 / sig))) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dcl[i] = (-Cov1[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    ddl[i] = (-Cov2[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    del[i] = (-Cov3[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dfl[i] = (-Cov4[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dgl[i] = (-Cov5[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dhl[i] = (-Cov6[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dil[i] = (-Cov7[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    djl[i] = (-Cov8[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dkl[i] = (-Cov9[i] * Cov10[i] * Math.Pow(time[i], (1 / sig)) * (Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig))) / (Math.Pow(sig, 2));
                                    dml[i] = ((censor[i] * Cov10[i]) / (Math.Pow(sig, 2))) - ((Cov10[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dm[i] = ((-censor[i] / sig) + ((censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 2))) - ((censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) + (((Math.Pow(time[i], (1 / sig)))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 2)));
                                    d2m[i] = (censor[i] / (Math.Pow(sig, 2))) - ((2 * censor[i] * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])) / (Math.Pow(sig, 3))) + ((2 * censor[i] * Math.Log(time[i])) / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * (2 * sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 4));
                                    dma[i] = (censor[i] / (Math.Pow(sig, 2))) - ((Math.Pow(time[i], (1 / sig))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmb[i] = ((censor[i] * SNP[i]) / (Math.Pow(sig, 2))) - ((SNP[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmc[i] = ((censor[i] * Cov1[i]) / (Math.Pow(sig, 2))) - ((Cov1[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmd[i] = ((censor[i] * Cov2[i]) / (Math.Pow(sig, 2))) - ((Cov2[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dme[i] = ((censor[i] * Cov3[i]) / (Math.Pow(sig, 2))) - ((Cov3[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmf[i] = ((censor[i] * Cov4[i]) / (Math.Pow(sig, 2))) - ((Cov4[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmg[i] = ((censor[i] * Cov5[i]) / (Math.Pow(sig, 2))) - ((Cov5[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmh[i] = ((censor[i] * Cov6[i]) / (Math.Pow(sig, 2))) - ((Cov6[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmi[i] = ((censor[i] * Cov7[i]) / (Math.Pow(sig, 2))) - ((Cov7[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmj[i] = ((censor[i] * Cov8[i]) / (Math.Pow(sig, 2))) - ((Cov8[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                    dmk[i] = ((censor[i] * Cov9[i]) / (Math.Pow(sig, 2))) - ((Cov9[i] * (Math.Pow(time[i], (1 / sig)))) * (sig + Math.Log(time[i]) - b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) * Math.Exp((-b0 - b1 * SNP[i] - b2 * Cov1[i] - b3 * Cov2[i] - b4 * Cov3[i] - b5 * Cov4[i] - b6 * Cov5[i] - b7 * Cov6[i] - b8 * Cov7[i] - b9 * Cov8[i] - b10 * Cov9[i] - b11 * Cov10[i]) / sig)) / (Math.Pow(sig, 3));
                                }
                            }
                            sumda = da.Sum();
                            sumdb = db.Sum();
                            sumdc = dc.Sum();
                            sumdd = dd.Sum();
                            sumde = de.Sum();
                            sumdf = df.Sum();
                            sumdg = dg.Sum();
                            sumdh = dh.Sum();
                            sumdi = di.Sum();
                            sumdj = dj.Sum();
                            sumdk = dk.Sum();
                            sumdl = dl.Sum();
                            sumdm = dm.Sum();

                            sumd2a = d2a.Sum();
                            sumd2b = d2b.Sum();
                            sumd2c = d2c.Sum();
                            sumd2d = d2d.Sum();
                            sumd2e = d2e.Sum();
                            sumd2f = d2f.Sum();
                            sumd2g = d2g.Sum();
                            sumd2h = d2h.Sum();
                            sumd2i = d2i.Sum();
                            sumd2j = d2j.Sum();
                            sumd2k = d2k.Sum();
                            sumd2l = d2l.Sum();
                            sumd2m = d2m.Sum();

                            sumdab = dab.Sum();
                            sumdac = dac.Sum();
                            sumdad = dad.Sum();
                            sumdae = dae.Sum();
                            sumdaf = daf.Sum();
                            sumdag = dag.Sum();
                            sumdah = dah.Sum();
                            sumdai = dai.Sum();
                            sumdaj = daj.Sum();
                            sumdak = dak.Sum();
                            sumdal = dal.Sum();

                            sumdbc = dbc.Sum();
                            sumdbd = dbd.Sum();
                            sumdbe = dbe.Sum();
                            sumdbf = dbf.Sum();
                            sumdbg = dbg.Sum();
                            sumdbh = dbh.Sum();
                            sumdbi = dbi.Sum();
                            sumdbj = dbj.Sum();
                            sumdbk = dbk.Sum();
                            sumdbl = dbl.Sum();

                            sumdcd = dcd.Sum();
                            sumdce = dce.Sum();
                            sumdcf = dcf.Sum();
                            sumdcg = dcg.Sum();
                            sumdch = dch.Sum();
                            sumdci = dci.Sum();
                            sumdcj = dcj.Sum();
                            sumdck = dck.Sum();
                            sumdcl = dcl.Sum();

                            sumdde = dde.Sum();
                            sumddf = ddf.Sum();
                            sumddg = ddg.Sum();
                            sumddh = ddh.Sum();
                            sumddi = ddi.Sum();
                            sumddj = ddj.Sum();
                            sumddk = ddk.Sum();
                            sumddl = ddl.Sum();

                            sumdef = def.Sum();
                            sumdeg = deg.Sum();
                            sumdeh = deh.Sum();
                            sumdei = dei.Sum();
                            sumdej = dej.Sum();
                            sumdek = dek.Sum();
                            sumdel = del.Sum();

                            sumdfg = dfg.Sum();
                            sumdfh = dfh.Sum();
                            sumdfi = dfi.Sum();
                            sumdfj = dfj.Sum();
                            sumdfk = dfk.Sum();
                            sumdfl = dfl.Sum();

                            sumdgh = dgh.Sum();
                            sumdgi = dgi.Sum();
                            sumdgj = dgj.Sum();
                            sumdgk = dgk.Sum();
                            sumdgl = dgl.Sum();

                            sumdhi = dhi.Sum();
                            sumdhj = dhj.Sum();
                            sumdhk = dhk.Sum();
                            sumdhl = dhl.Sum();

                            sumdij = dij.Sum();
                            sumdik = dik.Sum();
                            sumdil = dil.Sum();

                            sumdjk = djk.Sum();
                            sumdjl = djl.Sum();

                            sumdkl = dkl.Sum();

                            sumdma = dma.Sum();
                            sumdmb = dmb.Sum();
                            sumdmc = dmc.Sum();
                            sumdmd = dmd.Sum();
                            sumdme = dme.Sum();
                            sumdmf = dmf.Sum();
                            sumdmg = dmg.Sum();
                            sumdmh = dmh.Sum();
                            sumdmi = dmi.Sum();
                            sumdmj = dmj.Sum();
                            sumdmk = dmk.Sum();
                            sumdml = dml.Sum();

                            H1 = new double[13, 13] {{sumd2a,sumdab,sumdac,sumdad,sumdae,sumdaf,sumdag,sumdah,sumdai,sumdaj,sumdak,sumdal,sumdma},
                                                                {sumdab,sumd2b,sumdbc,sumdbd,sumdbe,sumdbf,sumdbg,sumdbh,sumdbi,sumdbj,sumdbk,sumdbl,sumdmb},
                                                                {sumdac,sumdbc,sumd2c,sumdcd,sumdce,sumdcf,sumdcg,sumdch,sumdci,sumdcj,sumdck,sumdcl,sumdmc},
                                                                {sumdad,sumdbd,sumdcd,sumd2d,sumdde,sumddf,sumddg,sumddh,sumddi,sumddj,sumddk,sumddl,sumdmd},
                                                                {sumdae,sumdbe,sumdce,sumdde,sumd2e,sumdef,sumdeg,sumdeh,sumdei,sumdej,sumdek,sumdel,sumdme},
                                                                {sumdaf,sumdbf,sumdcf,sumddf,sumdef,sumd2f,sumdfg,sumdfh,sumdfi,sumdfj,sumdfk,sumdfl,sumdmf},
                                                                {sumdag,sumdbg,sumdcg,sumddg,sumdeg,sumdfg,sumd2g,sumdgh,sumdgi,sumdgj,sumdgk,sumdgl,sumdmg},
                                                                {sumdah,sumdbh,sumdch,sumddh,sumdeh,sumdfh,sumdgh,sumd2h,sumdhi,sumdhj,sumdhk,sumdhl,sumdmh},
                                                                {sumdai,sumdbi,sumdci,sumddi,sumdei,sumdfi,sumdgi,sumdhi,sumd2i,sumdij,sumdik,sumdil,sumdmi},
                                                                {sumdaj,sumdbj,sumdcj,sumddj,sumdej,sumdfj,sumdgj,sumdhj,sumdij,sumd2j,sumdjk,sumdjl,sumdmj},
                                                                {sumdak,sumdbk,sumdck,sumddk,sumdek,sumdfk,sumdgk,sumdhk,sumdik,sumdjk,sumd2k,sumdkl,sumdmk},
                                                                {sumdal,sumdbl,sumdcl,sumddl,sumdel,sumdfl,sumdgl,sumdhl,sumdil,sumdjl,sumdkl,sumd2l,sumdml},
                                                                {sumdma,sumdmb,sumdmc,sumdmd,sumdme,sumdmf,sumdmg,sumdmh,sumdmi,sumdmj,sumdmk,sumdml,sumd2m}};



                            newArray2 = new double[covcount + intcount + 3, covcount + intcount + 3];
                            newArray2[0, 0] = H1[0, 0];
                            newArray2[0, 1] = H1[0, 1];
                            newArray2[1, 0] = H1[1, 0];
                            newArray2[1, 1] = H1[1, 1];

                            if (covcount + intcount + 3 >= 3)
                            {
                                newArray2[0, 2] = H1[0, 2];
                                newArray2[1, 2] = H1[1, 2];
                                newArray2[2, 2] = H1[2, 2];
                                newArray2[2, 0] = H1[2, 0];
                                newArray2[2, 1] = H1[2, 1];
                            }
                            if (covcount + intcount + 3 >= 4)
                            {
                                newArray2[0, 3] = H1[0, 3];
                                newArray2[1, 3] = H1[1, 3];
                                newArray2[2, 3] = H1[2, 3];
                                newArray2[3, 3] = H1[3, 3];
                                newArray2[3, 0] = H1[3, 0];
                                newArray2[3, 1] = H1[3, 1];
                                newArray2[3, 2] = H1[3, 2];
                            }
                            if (covcount + intcount + 3 >= 5)
                            {
                                newArray2[0, 4] = H1[0, 4];
                                newArray2[1, 4] = H1[1, 4];
                                newArray2[2, 4] = H1[2, 4];
                                newArray2[3, 4] = H1[3, 4];
                                newArray2[4, 4] = H1[4, 4];
                                newArray2[4, 0] = H1[4, 0];
                                newArray2[4, 1] = H1[4, 1];
                                newArray2[4, 2] = H1[4, 2];
                                newArray2[4, 3] = H1[4, 3];
                            }
                            if (covcount + intcount + 3 >= 6)
                            {
                                newArray2[0, 5] = H1[0, 5];
                                newArray2[1, 5] = H1[1, 5];
                                newArray2[2, 5] = H1[2, 5];
                                newArray2[3, 5] = H1[3, 5];
                                newArray2[4, 5] = H1[4, 5];
                                newArray2[5, 5] = H1[5, 5];
                                newArray2[5, 0] = H1[5, 0];
                                newArray2[5, 1] = H1[5, 1];
                                newArray2[5, 2] = H1[5, 2];
                                newArray2[5, 3] = H1[5, 3];
                                newArray2[5, 4] = H1[5, 4];
                            }
                            if (covcount + intcount + 3 >= 7)
                            {
                                newArray2[0, 6] = H1[0, 6];
                                newArray2[1, 6] = H1[1, 6];
                                newArray2[2, 6] = H1[2, 6];
                                newArray2[3, 6] = H1[3, 6];
                                newArray2[4, 6] = H1[4, 6];
                                newArray2[5, 6] = H1[5, 6];
                                newArray2[6, 6] = H1[6, 6];
                                newArray2[6, 0] = H1[6, 0];
                                newArray2[6, 1] = H1[6, 1];
                                newArray2[6, 2] = H1[6, 2];
                                newArray2[6, 3] = H1[6, 3];
                                newArray2[6, 4] = H1[6, 4];
                                newArray2[6, 5] = H1[6, 5];
                            }
                            if (covcount + intcount + 3 >= 8)
                            {
                                newArray2[0, 7] = H1[0, 7];
                                newArray2[1, 7] = H1[1, 7];
                                newArray2[2, 7] = H1[2, 7];
                                newArray2[3, 7] = H1[3, 7];
                                newArray2[4, 7] = H1[4, 7];
                                newArray2[5, 7] = H1[5, 7];
                                newArray2[6, 7] = H1[6, 7];
                                newArray2[7, 7] = H1[7, 7];
                                newArray2[7, 0] = H1[7, 0];
                                newArray2[7, 1] = H1[7, 1];
                                newArray2[7, 2] = H1[7, 2];
                                newArray2[7, 3] = H1[7, 3];
                                newArray2[7, 4] = H1[7, 4];
                                newArray2[7, 5] = H1[7, 5];
                                newArray2[7, 6] = H1[7, 6];
                            }
                            if (covcount + intcount + 3 >= 9)
                            {
                                newArray2[0, 8] = H1[0, 8];
                                newArray2[1, 8] = H1[1, 8];
                                newArray2[2, 8] = H1[2, 8];
                                newArray2[3, 8] = H1[3, 8];
                                newArray2[4, 8] = H1[4, 8];
                                newArray2[5, 8] = H1[5, 8];
                                newArray2[6, 8] = H1[6, 8];
                                newArray2[7, 8] = H1[7, 8];
                                newArray2[8, 8] = H1[8, 8];
                                newArray2[8, 0] = H1[8, 0];
                                newArray2[8, 1] = H1[8, 1];
                                newArray2[8, 2] = H1[8, 2];
                                newArray2[8, 3] = H1[8, 3];
                                newArray2[8, 4] = H1[8, 4];
                                newArray2[8, 5] = H1[8, 5];
                                newArray2[8, 6] = H1[8, 6];
                                newArray2[8, 7] = H1[8, 7];
                            }
                            if (covcount + intcount + 3 >= 10)
                            {
                                newArray2[0, 9] = H1[0, 9];
                                newArray2[1, 9] = H1[1, 9];
                                newArray2[2, 9] = H1[2, 9];
                                newArray2[3, 9] = H1[3, 9];
                                newArray2[4, 9] = H1[4, 9];
                                newArray2[5, 9] = H1[5, 9];
                                newArray2[6, 9] = H1[6, 9];
                                newArray2[7, 9] = H1[7, 9];
                                newArray2[8, 9] = H1[8, 9];
                                newArray2[9, 9] = H1[9, 9];
                                newArray2[9, 0] = H1[9, 0];
                                newArray2[9, 1] = H1[9, 1];
                                newArray2[9, 2] = H1[9, 2];
                                newArray2[9, 3] = H1[9, 3];
                                newArray2[9, 4] = H1[9, 4];
                                newArray2[9, 5] = H1[9, 5];
                                newArray2[9, 6] = H1[9, 6];
                                newArray2[9, 7] = H1[9, 7];
                                newArray2[9, 8] = H1[9, 8];
                            }
                            if (covcount + intcount + 3 >= 11)
                            {
                                newArray2[0, 10] = H1[0, 10];
                                newArray2[1, 10] = H1[1, 10];
                                newArray2[2, 10] = H1[2, 10];
                                newArray2[3, 10] = H1[3, 10];
                                newArray2[4, 10] = H1[4, 10];
                                newArray2[5, 10] = H1[5, 10];
                                newArray2[6, 10] = H1[6, 10];
                                newArray2[7, 10] = H1[7, 10];
                                newArray2[8, 10] = H1[8, 10];
                                newArray2[9, 10] = H1[9, 10];
                                newArray2[10, 10] = H1[10, 10];
                                newArray2[10, 0] = H1[10, 0];
                                newArray2[10, 1] = H1[10, 1];
                                newArray2[10, 2] = H1[10, 2];
                                newArray2[10, 3] = H1[10, 3];
                                newArray2[10, 4] = H1[10, 4];
                                newArray2[10, 5] = H1[10, 5];
                                newArray2[10, 6] = H1[10, 6];
                                newArray2[10, 7] = H1[10, 7];
                                newArray2[10, 8] = H1[10, 8];
                                newArray2[10, 9] = H1[10, 9];
                            }
                            if (covcount + intcount + 3 >= 12)
                            {
                                newArray2[0, 11] = H1[0, 11];
                                newArray2[1, 11] = H1[1, 11];
                                newArray2[2, 11] = H1[2, 11];
                                newArray2[3, 11] = H1[3, 11];
                                newArray2[4, 11] = H1[4, 11];
                                newArray2[5, 11] = H1[5, 11];
                                newArray2[6, 11] = H1[6, 11];
                                newArray2[7, 11] = H1[7, 11];
                                newArray2[8, 11] = H1[8, 11];
                                newArray2[9, 11] = H1[9, 11];
                                newArray2[10, 11] = H1[10, 11];
                                newArray2[11, 0] = H1[11, 0];
                                newArray2[11, 1] = H1[11, 1];
                                newArray2[11, 2] = H1[11, 2];
                                newArray2[11, 3] = H1[11, 3];
                                newArray2[11, 4] = H1[11, 4];
                                newArray2[11, 5] = H1[11, 5];
                                newArray2[11, 6] = H1[11, 6];
                                newArray2[11, 7] = H1[11, 7];
                                newArray2[11, 8] = H1[11, 8];
                                newArray2[11, 9] = H1[11, 9];
                                newArray2[11, 10] = H1[11, 10];
                                newArray2[11, 11] = H1[11, 11];
                            }
                            if (covcount + intcount + 3 >= 13)
                            {
                                newArray2[0, 12] = H1[0, 12];
                                newArray2[1, 12] = H1[1, 12];
                                newArray2[2, 12] = H1[2, 12];
                                newArray2[3, 12] = H1[3, 12];
                                newArray2[4, 12] = H1[4, 12];
                                newArray2[5, 12] = H1[5, 12];
                                newArray2[6, 12] = H1[6, 12];
                                newArray2[7, 12] = H1[7, 12];
                                newArray2[8, 12] = H1[8, 12];
                                newArray2[9, 12] = H1[9, 12];
                                newArray2[10, 12] = H1[10, 12];
                                newArray2[11, 12] = H1[11, 12];
                                newArray2[12, 0] = H1[12, 0];
                                newArray2[12, 1] = H1[12, 1];
                                newArray2[12, 2] = H1[12, 2];
                                newArray2[12, 3] = H1[12, 3];
                                newArray2[12, 4] = H1[12, 4];
                                newArray2[12, 5] = H1[12, 5];
                                newArray2[12, 6] = H1[12, 6];
                                newArray2[12, 7] = H1[12, 7];
                                newArray2[12, 8] = H1[12, 8];
                                newArray2[12, 9] = H1[12, 9];
                                newArray2[12, 10] = H1[12, 10];
                                newArray2[12, 11] = H1[12, 11];
                                newArray2[12, 12] = H1[12, 12];
                            }

                            H1 = newArray2;

                            IH = new double[covcount + intcount + 3, covcount + intcount + 3];
                            try
                            {
                                IH = H1.Inverse();
                            }
                            catch (Exception next)
                            {

                                valid = 1;
                                goto Try;
                            }
                            valid = 0;


                            Z = new double[,]{{sumda},
                                                       {sumdb},
                                                       {sumdc},
                                                       {sumdd},
                                                       {sumde},
                                                       {sumdf},
                                                       {sumdg},
                                                       {sumdh},
                                                       {sumdi},
                                                       {sumdj},
                                                       {sumdk},
                                                       {sumdl},
                                                       {sumdm}};

                            newArray4 = new double[covcount + intcount + 3, 1];
                            Array.Copy(Z, newArray4, newArray4.Length);
                            Z = newArray4;

                            // newton raphson

                            NR = ini.Subtract(IH.Multiply(Z));

                            if (NR[0, 0].Equals(Double.NaN) || NR[1, 0].Equals(Double.NaN))
                            {
                                valid = 1;
                                goto Try;
                             // break;
                            }
                      /*      var weiblike = new double[dos.Rows.Count];
                            var weiblike2 = new double[dos.Rows.Count];
                            for (int i = 0; i < dos.Rows.Count; ++i)
                            {
                                weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + 0 * Cov6[i] + 0 * Cov7[i] + 0 * Cov8[i] + 0 * Cov9[i] + 0 * Cov10[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + 0 * Cov6[i] + 0 * Cov7[i] + 0 * Cov8[i] + 0 * Cov9[i] + 0 * Cov10[i])))));
                                weiblike2[i] = ((censor[i]) * (-Math.Log(sig) + (((1 / sig) - 1) * Math.Log(time[i])) - ((1 / sig) * (b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5* Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i]))) - (time[i] / (sig * (Math.Exp(b0 + b1 * SNP[i] + b2 * Cov1[i] + b3 * Cov2[i] + b4 * Cov3[i] + b5 * Cov4[i] + b6 * Cov5[i] + b7 * Cov6[i] + b8 * Cov7[i] + b9 * Cov8[i] + b10 * Cov9[i] + b11 * Cov10[i])))));
                            }
                            double loglikealt = weiblike.Sum();
                            double loglikestart = weiblike2.Sum();
                            
                            */

                            if (Math.Abs(NR[1, 0] - b1) < 0.000001 && Math.Abs(NR[covcount + intcount + 2, 0] - sig) < 0.001)
                            // Math.Abs(loglikealt - loglikestart) < 0.1
                            { break; }
                            else
                            {
                                continue;
                            }
                        }

                        if (NR[0, 0].Equals(Double.NaN) || NR[1, 0].Equals(Double.NaN))
                        {
                            Console.WriteLine("NaN's produced, " + dos.Columns[0] + " analysis did not converge");
                            goto SKIP2;
                        }

                        //variance

                        V = new double[13, 13] {{-sumd2a,-sumdab,-sumdac,-sumdad,-sumdae,-sumdaf,-sumdag,-sumdah,-sumdai,-sumdaj,-sumdak,-sumdal,-sumdma},
                                                                    {-sumdab,-sumd2b,-sumdbc,-sumdbd,-sumdbe,-sumdbf,-sumdbg,-sumdbh,-sumdbi,-sumdbj,-sumdbk,-sumdbl,-sumdmb},
                                                                    {-sumdac,-sumdbc,-sumd2c,-sumdcd,-sumdce,-sumdcf,-sumdcg,-sumdch,-sumdci,-sumdcj,-sumdck,-sumdcl,-sumdmc},
                                                                    {-sumdad,-sumdbd,-sumdcd,-sumd2d,-sumdde,-sumddf,-sumddg,-sumddh,-sumddi,-sumddj,-sumddk,-sumddl,-sumdmd},
                                                                    {-sumdae,-sumdbe,-sumdce,-sumdde,-sumd2e,-sumdef,-sumdeg,-sumdeh,-sumdei,-sumdej,-sumdek,-sumdel,-sumdme},
                                                                    {-sumdaf,-sumdbf,-sumdcf,-sumddf,-sumdef,-sumd2f,-sumdfg,-sumdfh,-sumdfi,-sumdfj,-sumdfk,-sumdfl,-sumdmf},
                                                                    {-sumdag,-sumdbg,-sumdcg,-sumddg,-sumdeg,-sumdfg,-sumd2g,-sumdgh,-sumdgi,-sumdgj,-sumdgk,-sumdgl,-sumdmg},
                                                                    {-sumdah,-sumdbh,-sumdch,-sumddh,-sumdeh,-sumdfh,-sumdgh,-sumd2h,-sumdhi,-sumdhj,-sumdhk,-sumdhl,-sumdmh},
                                                                    {-sumdai,-sumdbi,-sumdci,-sumddi,-sumdei,-sumdfi,-sumdgi,-sumdhi,-sumd2i,-sumdij,-sumdik,-sumdil,-sumdmi},
                                                                    {-sumdaj,-sumdbj,-sumdcj,-sumddj,-sumdej,-sumdfj,-sumdgj,-sumdhj,-sumdij,-sumd2j,-sumdjk,-sumdjl,-sumdmj},
                                                                    {-sumdak,-sumdbk,-sumdck,-sumddk,-sumdek,-sumdfk,-sumdgk,-sumdhk,-sumdik,-sumdjk,-sumd2k,-sumdkl,-sumdmk},
                                                                    {-sumdal,-sumdbl,-sumdcl,-sumddl,-sumdel,-sumdfl,-sumdgl,-sumdhl,-sumdil,-sumdjl,-sumdkl,-sumd2l,-sumdml},
                                                                    {-sumdma,-sumdmb,-sumdmc,-sumdmd,-sumdme,-sumdmf,-sumdmg,-sumdmh,-sumdmi,-sumdmj,-sumdmk,-sumdml,-sumd2m}};

                        newArray5 = new double[covcount + intcount + 3, covcount + intcount + 3];
                        newArray5[0, 0] = V[0, 0];
                        newArray5[0, 1] = V[0, 1];
                        newArray5[1, 0] = V[1, 0];
                        newArray5[1, 1] = V[1, 1];

                        if (covcount + intcount + 3 >= 3)
                        {
                            newArray5[0, 2] = V[0, 2];
                            newArray5[1, 2] = V[1, 2];
                            newArray5[2, 2] = V[2, 2];
                            newArray5[2, 0] = V[2, 0];
                            newArray5[2, 1] = V[2, 1];
                        }
                        if (covcount + intcount + 3 >= 4)
                        {
                            newArray5[0, 3] = V[0, 3];
                            newArray5[1, 3] = V[1, 3];
                            newArray5[2, 3] = V[2, 3];
                            newArray5[3, 3] = V[3, 3];
                            newArray5[3, 0] = V[3, 0];
                            newArray5[3, 1] = V[3, 1];
                            newArray5[3, 2] = V[3, 2];
                        }
                        if (covcount + intcount + 3 >= 5)
                        {
                            newArray5[0, 4] = V[0, 4];
                            newArray5[1, 4] = V[1, 4];
                            newArray5[2, 4] = V[2, 4];
                            newArray5[3, 4] = V[3, 4];
                            newArray5[4, 4] = V[4, 4];
                            newArray5[4, 0] = V[4, 0];
                            newArray5[4, 1] = V[4, 1];
                            newArray5[4, 2] = V[4, 2];
                            newArray5[4, 3] = V[4, 3];
                        }
                        if (covcount + intcount + 3 >= 6)
                        {
                            newArray5[0, 5] = V[0, 5];
                            newArray5[1, 5] = V[1, 5];
                            newArray5[2, 5] = V[2, 5];
                            newArray5[3, 5] = V[3, 5];
                            newArray5[4, 5] = V[4, 5];
                            newArray5[5, 5] = V[5, 5];
                            newArray5[5, 0] = V[5, 0];
                            newArray5[5, 1] = V[5, 1];
                            newArray5[5, 2] = V[5, 2];
                            newArray5[5, 3] = V[5, 3];
                            newArray5[5, 4] = V[5, 4];
                        }
                        if (covcount + intcount + 3 >= 7)
                        {
                            newArray5[0, 6] = V[0, 6];
                            newArray5[1, 6] = V[1, 6];
                            newArray5[2, 6] = V[2, 6];
                            newArray5[3, 6] = V[3, 6];
                            newArray5[4, 6] = V[4, 6];
                            newArray5[5, 6] = V[5, 6];
                            newArray5[6, 6] = V[6, 6];
                            newArray5[6, 0] = V[6, 0];
                            newArray5[6, 1] = V[6, 1];
                            newArray5[6, 2] = V[6, 2];
                            newArray5[6, 3] = V[6, 3];
                            newArray5[6, 4] = V[6, 4];
                            newArray5[6, 5] = V[6, 5];
                        }
                        if (covcount + intcount + 3 >= 8)
                        {
                            newArray5[0, 7] = V[0, 7];
                            newArray5[1, 7] = V[1, 7];
                            newArray5[2, 7] = V[2, 7];
                            newArray5[3, 7] = V[3, 7];
                            newArray5[4, 7] = V[4, 7];
                            newArray5[5, 7] = V[5, 7];
                            newArray5[6, 7] = V[6, 7];
                            newArray5[7, 7] = V[7, 7];
                            newArray5[7, 0] = V[7, 0];
                            newArray5[7, 1] = V[7, 1];
                            newArray5[7, 2] = V[7, 2];
                            newArray5[7, 3] = V[7, 3];
                            newArray5[7, 4] = V[7, 4];
                            newArray5[7, 5] = V[7, 5];
                            newArray5[7, 6] = V[7, 6];
                        }
                        if (covcount + intcount + 3 >= 9)
                        {
                            newArray5[0, 8] = V[0, 8];
                            newArray5[1, 8] = V[1, 8];
                            newArray5[2, 8] = V[2, 8];
                            newArray5[3, 8] = V[3, 8];
                            newArray5[4, 8] = V[4, 8];
                            newArray5[5, 8] = V[5, 8];
                            newArray5[6, 8] = V[6, 8];
                            newArray5[7, 8] = V[7, 8];
                            newArray5[8, 8] = V[8, 8];
                            newArray5[8, 0] = V[8, 0];
                            newArray5[8, 1] = V[8, 1];
                            newArray5[8, 2] = V[8, 2];
                            newArray5[8, 3] = V[8, 3];
                            newArray5[8, 4] = V[8, 4];
                            newArray5[8, 5] = V[8, 5];
                            newArray5[8, 6] = V[8, 6];
                            newArray5[8, 7] = V[8, 7];
                        }
                        if (covcount + intcount + 3 >= 10)
                        {
                            newArray5[0, 9] = V[0, 9];
                            newArray5[1, 9] = V[1, 9];
                            newArray5[2, 9] = V[2, 9];
                            newArray5[3, 9] = V[3, 9];
                            newArray5[4, 9] = V[4, 9];
                            newArray5[5, 9] = V[5, 9];
                            newArray5[6, 9] = V[6, 9];
                            newArray5[7, 9] = V[7, 9];
                            newArray5[8, 9] = V[8, 9];
                            newArray5[9, 9] = V[9, 9];
                            newArray5[9, 0] = V[9, 0];
                            newArray5[9, 1] = V[9, 1];
                            newArray5[9, 2] = V[9, 2];
                            newArray5[9, 3] = V[9, 3];
                            newArray5[9, 4] = V[9, 4];
                            newArray5[9, 5] = V[9, 5];
                            newArray5[9, 6] = V[9, 6];
                            newArray5[9, 7] = V[9, 7];
                            newArray5[9, 8] = V[9, 8];
                        }
                        if (covcount + intcount + 3 >= 11)
                        {
                            newArray5[0, 10] = V[0, 10];
                            newArray5[1, 10] = V[1, 10];
                            newArray5[2, 10] = V[2, 10];
                            newArray5[3, 10] = V[3, 10];
                            newArray5[4, 10] = V[4, 10];
                            newArray5[5, 10] = V[5, 10];
                            newArray5[6, 10] = V[6, 10];
                            newArray5[7, 10] = V[7, 10];
                            newArray5[8, 10] = V[8, 10];
                            newArray5[9, 10] = V[9, 10];
                            newArray5[10, 10] = V[10, 10];
                            newArray5[10, 0] = V[10, 0];
                            newArray5[10, 1] = V[10, 1];
                            newArray5[10, 2] = V[10, 2];
                            newArray5[10, 3] = V[10, 3];
                            newArray5[10, 4] = V[10, 4];
                            newArray5[10, 5] = V[10, 5];
                            newArray5[10, 6] = V[10, 6];
                            newArray5[10, 7] = V[10, 7];
                            newArray5[10, 8] = V[10, 8];
                            newArray5[10, 9] = V[10, 9];
                        }
                        if (covcount + intcount + 3 >= 12)
                        {
                            newArray5[0, 11] = V[0, 11];
                            newArray5[1, 11] = V[1, 11];
                            newArray5[2, 11] = V[2, 11];
                            newArray5[3, 11] = V[3, 11];
                            newArray5[4, 11] = V[4, 11];
                            newArray5[5, 11] = V[5, 11];
                            newArray5[6, 11] = V[6, 11];
                            newArray5[7, 11] = V[7, 11];
                            newArray5[8, 11] = V[8, 11];
                            newArray5[9, 11] = V[9, 11];
                            newArray5[10, 11] = V[10, 11];
                            newArray5[11, 0] = V[11, 0];
                            newArray5[11, 1] = V[11, 1];
                            newArray5[11, 2] = V[11, 2];
                            newArray5[11, 3] = V[11, 3];
                            newArray5[11, 4] = V[11, 4];
                            newArray5[11, 5] = V[11, 5];
                            newArray5[11, 6] = V[11, 6];
                            newArray5[11, 7] = V[11, 7];
                            newArray5[11, 8] = V[11, 8];
                            newArray5[11, 9] = V[11, 9];
                            newArray5[11, 10] = V[11, 10];
                            newArray5[11, 11] = V[11, 11];
                        }
                        if (covcount + intcount + 3 >= 13)
                        {
                            newArray5[0, 12] = V[0, 12];
                            newArray5[1, 12] = V[1, 12];
                            newArray5[2, 12] = V[2, 12];
                            newArray5[3, 12] = V[3, 12];
                            newArray5[4, 12] = V[4, 12];
                            newArray5[5, 12] = V[5, 12];
                            newArray5[6, 12] = V[6, 12];
                            newArray5[7, 12] = V[7, 12];
                            newArray5[8, 12] = V[8, 12];
                            newArray5[9, 12] = V[9, 12];
                            newArray5[10, 12] = V[10, 12];
                            newArray5[11, 12] = V[11, 12];
                            newArray5[12, 0] = V[12, 0];
                            newArray5[12, 1] = V[12, 1];
                            newArray5[12, 2] = V[12, 2];
                            newArray5[12, 3] = V[12, 3];
                            newArray5[12, 4] = V[12, 4];
                            newArray5[12, 5] = V[12, 5];
                            newArray5[12, 6] = V[12, 6];
                            newArray5[12, 7] = V[12, 7];
                            newArray5[12, 8] = V[12, 8];
                            newArray5[12, 9] = V[12, 9];
                            newArray5[12, 10] = V[12, 10];
                            newArray5[12, 11] = V[12, 11];
                            newArray5[12, 12] = V[12, 12];
                        }

                        V = newArray5;

                        IV = V.Inverse();

                        SE = (IV).Sqrt();

                        // z score

                        z0 = NR[0, 0] / SE[0, 0];
                        z1 = NR[1, 0] / SE[1, 1];

                        Func<double, double> g = (x) => (1 / Math.Sqrt(2 * Math.PI)) * Math.Exp(-(x * x) / 2);

                        if (covcount + intcount + 3 == 3)
                        {
                            z2 = NR[2, 0] / SE[2, 2];

                            zscore = new double[,] { { z0 }, { z1 }, { z2 } };

                            // p- values
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval } };

                            //waldtest
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);

                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue } };

                        }


                        if (covcount + intcount + 3 == 4)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];

                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }

                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval } };

                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue } };
                           
                        }

                        if (covcount + intcount + 3 == 5)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }

                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval } };

                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue } };

                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);

                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 6)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }

                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue } };
                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);
                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 7)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            z6 = NR[6, 0] / SE[6, 6];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 }, { z6 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }
                            if (z6 > 0)
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z6);
                            }
                            else
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, z6, Double.PositiveInfinity);
                            }

                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);
                            b6pval = 2 * (1 - iagk6);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval }, { b6pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            w6 = new WaldTest(z6);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue }, { w6.PValue } };
                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);
                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 8)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            z6 = NR[6, 0] / SE[6, 6];
                            z7 = NR[7, 0] / SE[7, 7];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 }, { z6 }, { z7 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }
                            if (z6 > 0)
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z6);
                            }
                            else
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, z6, Double.PositiveInfinity);
                            }
                            if (z7 > 0)
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z7);
                            }
                            else
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, z7, Double.PositiveInfinity);
                            }

                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);
                            b6pval = 2 * (1 - iagk6);
                            b7pval = 2 * (1 - iagk7);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval }, { b6pval }, { b7pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            w6 = new WaldTest(z6);
                            w7 = new WaldTest(z7);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue }, { w6.PValue }, { w7.PValue } };
                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i]+ NR[5, 0] * Cov4[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);
                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 9)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            z6 = NR[6, 0] / SE[6, 6];
                            z7 = NR[7, 0] / SE[7, 7];
                            z8 = NR[8, 0] / SE[8, 8];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 }, { z6 }, { z7 }, { z8 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }
                            if (z6 > 0)
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z6);
                            }
                            else
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, z6, Double.PositiveInfinity);
                            }
                            if (z7 > 0)
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z7);
                            }
                            else
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, z7, Double.PositiveInfinity);
                            }
                            if (z8 > 0)
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z8);
                            }
                            else
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, z8, Double.PositiveInfinity);
                            }


                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);
                            b6pval = 2 * (1 - iagk6);
                            b7pval = 2 * (1 - iagk7);
                            b8pval = 2 * (1 - iagk8);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval }, { b6pval }, { b7pval }, { b8pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            w6 = new WaldTest(z6);
                            w7 = new WaldTest(z7);
                            w8 = new WaldTest(z8);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue }, { w6.PValue }, { w7.PValue }, { w8.PValue } };
                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);
                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 10)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            z6 = NR[6, 0] / SE[6, 6];
                            z7 = NR[7, 0] / SE[7, 7];
                            z8 = NR[8, 0] / SE[8, 8];
                            z9 = NR[9, 0] / SE[9, 9];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 }, { z6 }, { z7 }, { z8 }, { z9 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }
                            if (z6 > 0)
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z6);
                            }
                            else
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, z6, Double.PositiveInfinity);
                            }
                            if (z7 > 0)
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z7);
                            }
                            else
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, z7, Double.PositiveInfinity);
                            }
                            if (z8 > 0)
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z8);
                            }
                            else
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, z8, Double.PositiveInfinity);
                            }
                            if (z9 > 0)
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z9);
                            }
                            else
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, z9, Double.PositiveInfinity);
                            }

                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);
                            b6pval = 2 * (1 - iagk6);
                            b7pval = 2 * (1 - iagk7);
                            b8pval = 2 * (1 - iagk8);
                            b9pval = 2 * (1 - iagk9);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval }, { b6pval }, { b7pval }, { b8pval }, { b9pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            w6 = new WaldTest(z6);
                            w7 = new WaldTest(z7);
                            w8 = new WaldTest(z8);
                            w9 = new WaldTest(z9);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue }, { w6.PValue }, { w7.PValue }, { w8.PValue }, { w9.PValue } };

                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);
                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 11)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            z6 = NR[6, 0] / SE[6, 6];
                            z7 = NR[7, 0] / SE[7, 7];
                            z8 = NR[8, 0] / SE[8, 8];
                            z9 = NR[9, 0] / SE[9, 9];
                            z10 = NR[10, 0] / SE[10, 10];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 }, { z6 }, { z7 }, { z8 }, { z9 }, { z10 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }
                            if (z6 > 0)
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z6);
                            }
                            else
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, z6, Double.PositiveInfinity);
                            }
                            if (z7 > 0)
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z7);
                            }
                            else
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, z7, Double.PositiveInfinity);
                            }
                            if (z8 > 0)
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z8);
                            }
                            else
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, z8, Double.PositiveInfinity);
                            }
                            if (z9 > 0)
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z9);
                            }
                            else
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, z9, Double.PositiveInfinity);
                            }
                            if (z10 > 0)
                            {
                                iagk10 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z10);
                            }
                            else
                            {
                                iagk10 = InfiniteAdaptiveGaussKronrod.Integrate(g, z10, Double.PositiveInfinity);
                            }


                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);
                            b6pval = 2 * (1 - iagk6);
                            b7pval = 2 * (1 - iagk7);
                            b8pval = 2 * (1 - iagk8);
                            b9pval = 2 * (1 - iagk9);
                            b10pval = 2 * (1 - iagk10);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval }, { b6pval }, { b7pval }, { b8pval }, { b9pval }, { b10pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            w6 = new WaldTest(z6);
                            w7 = new WaldTest(z7);
                            w8 = new WaldTest(z8);
                            w9 = new WaldTest(z9);
                            w10 = new WaldTest(z10);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue }, { w6.PValue }, { w7.PValue }, { w8.PValue }, { w9.PValue }, { w10.PValue } };
                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i] + NR[9, 0] * Cov8[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i] + NR[9, 0] * Cov8[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);
                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 12)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            z6 = NR[6, 0] / SE[6, 6];
                            z7 = NR[7, 0] / SE[7, 7];
                            z8 = NR[8, 0] / SE[8, 8];
                            z9 = NR[9, 0] / SE[9, 9];
                            z10 = NR[10, 0] / SE[10, 10];
                            z11 = NR[11, 0] / SE[11, 11];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 }, { z6 }, { z7 }, { z8 }, { z9 }, { z10 }, { z11 } };
                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }
                            if (z6 > 0)
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z6);
                            }
                            else
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, z6, Double.PositiveInfinity);
                            }
                            if (z7 > 0)
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z7);
                            }
                            else
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, z7, Double.PositiveInfinity);
                            }
                            if (z8 > 0)
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z8);
                            }
                            else
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, z8, Double.PositiveInfinity);
                            }
                            if (z9 > 0)
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z9);
                            }
                            else
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, z9, Double.PositiveInfinity);
                            }
                            if (z10 > 0)
                            {
                                iagk10 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z10);
                            }
                            else
                            {
                                iagk10 = InfiniteAdaptiveGaussKronrod.Integrate(g, z10, Double.PositiveInfinity);
                            }
                            if (z11 > 0)
                            {
                                iagk11 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z11);
                            }
                            else
                            {
                                iagk11 = InfiniteAdaptiveGaussKronrod.Integrate(g, z11, Double.PositiveInfinity);
                            }

                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);
                            b6pval = 2 * (1 - iagk6);
                            b7pval = 2 * (1 - iagk7);
                            b8pval = 2 * (1 - iagk8);
                            b9pval = 2 * (1 - iagk9);
                            b10pval = 2 * (1 - iagk10);
                            b11pval = 2 * (1 - iagk11);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval }, { b6pval }, { b7pval }, { b8pval }, { b9pval }, { b10pval }, { b11pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            w6 = new WaldTest(z6);
                            w7 = new WaldTest(z7);
                            w8 = new WaldTest(z8);
                            w9 = new WaldTest(z9);
                            w10 = new WaldTest(z10);
                            w11 = new WaldTest(z11);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue }, { w6.PValue }, { w7.PValue }, { w8.PValue }, { w9.PValue }, { w10.PValue }, { w11.PValue } };
                            if (covariates != null & int1 != null)
                            {
                                weiblike = new double[dos.Rows.Count];
                                weiblike2 = new double[dos.Rows.Count];
                                for (int i = 0; i < dos.Rows.Count; ++i)
                                {
                                    weiblike[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i] + NR[9, 0] * Cov8[i] + NR[10, 0] * Cov9[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[1, 0] * SNP[i] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i] + NR[9, 0] * Cov8[i] + NR[10, 0] * Cov9[i])))));
                                    weiblike2[i] = ((censor[i]) * (-Math.Log(NR[covcount + intcount + 2, 0]) + (((1 / NR[covcount + intcount + 2, 0]) - 1) * Math.Log(time[i])) - ((1 / NR[covcount + intcount + 2, 0]) * (NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i] + NR[9, 0] * Cov8[i]))) - (time[i] / (NR[covcount + intcount + 2, 0] * (Math.Exp(NR[0, 0] + NR[2, 0] * Cov1[i] + NR[3, 0] * Cov2[i] + NR[4, 0] * Cov3[i] + NR[5, 0] * Cov4[i] + NR[6, 0] * Cov5[i] + NR[7, 0] * Cov6[i] + NR[8, 0] * Cov7[i] + NR[9, 0] * Cov8[i])))));
                                }
                                loglikealt = weiblike.Sum();
                                loglikenull = weiblike2.Sum();
                                likeratiotest = 2 * (loglikealt - loglikenull);
                                chi = new ChiSquareTest(likeratiotest, 2);
                                jointint = chi.PValue;
                            }
                        }
                        if (covcount + intcount + 3 == 13)
                        {
                            z2 = NR[2, 0] / SE[2, 2];
                            z3 = NR[3, 0] / SE[3, 3];
                            z4 = NR[4, 0] / SE[4, 4];
                            z5 = NR[5, 0] / SE[5, 5];
                            z6 = NR[6, 0] / SE[6, 6];
                            z7 = NR[7, 0] / SE[7, 7];
                            z8 = NR[8, 0] / SE[8, 8];
                            z9 = NR[9, 0] / SE[9, 9];
                            z10 = NR[10, 0] / SE[10, 10];
                            z11 = NR[11, 0] / SE[11, 11];
                            z12 = NR[12, 0] / SE[12, 12];
                            zscore = new double[,] { { z0 }, { z1 }, { z2 }, { z3 }, { z4 }, { z5 }, { z6 }, { z7 }, { z8 }, { z9 }, { z10 }, { z11 }, { z12 } };

                            if (z0 > 0)
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z0);
                            }
                            else
                            {
                                iagk0 = InfiniteAdaptiveGaussKronrod.Integrate(g, z0, Double.PositiveInfinity);
                            }
                            if (z1 > 0)
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z1);
                            }
                            else
                            {
                                iagk = InfiniteAdaptiveGaussKronrod.Integrate(g, z1, Double.PositiveInfinity);
                            }
                            if (z2 > 0)
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z2);
                            }
                            else
                            {
                                iagk2 = InfiniteAdaptiveGaussKronrod.Integrate(g, z2, Double.PositiveInfinity);
                            }
                            if (z3 > 0)
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z3);
                            }
                            else
                            {
                                iagk3 = InfiniteAdaptiveGaussKronrod.Integrate(g, z3, Double.PositiveInfinity);
                            }
                            if (z4 > 0)
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z4);
                            }
                            else
                            {
                                iagk4 = InfiniteAdaptiveGaussKronrod.Integrate(g, z4, Double.PositiveInfinity);
                            }
                            if (z5 > 0)
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z5);
                            }
                            else
                            {
                                iagk5 = InfiniteAdaptiveGaussKronrod.Integrate(g, z5, Double.PositiveInfinity);
                            }
                            if (z6 > 0)
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z6);
                            }
                            else
                            {
                                iagk6 = InfiniteAdaptiveGaussKronrod.Integrate(g, z6, Double.PositiveInfinity);
                            }
                            if (z7 > 0)
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z7);
                            }
                            else
                            {
                                iagk7 = InfiniteAdaptiveGaussKronrod.Integrate(g, z7, Double.PositiveInfinity);
                            }
                            if (z8 > 0)
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z8);
                            }
                            else
                            {
                                iagk8 = InfiniteAdaptiveGaussKronrod.Integrate(g, z8, Double.PositiveInfinity);
                            }
                            if (z9 > 0)
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z9);
                            }
                            else
                            {
                                iagk9 = InfiniteAdaptiveGaussKronrod.Integrate(g, z9, Double.PositiveInfinity);
                            }
                            if (z10 > 0)
                            {
                                iagk10 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z10);
                            }
                            else
                            {
                                iagk10 = InfiniteAdaptiveGaussKronrod.Integrate(g, z10, Double.PositiveInfinity);
                            }
                            if (z11 > 0)
                            {
                                iagk11 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z11);
                            }
                            else
                            {
                                iagk11 = InfiniteAdaptiveGaussKronrod.Integrate(g, z11, Double.PositiveInfinity);
                            }
                            if (z12 > 0)
                            {
                                iagk12 = InfiniteAdaptiveGaussKronrod.Integrate(g, Double.NegativeInfinity, z12);
                            }
                            else
                            {
                                iagk12 = InfiniteAdaptiveGaussKronrod.Integrate(g, z12, Double.PositiveInfinity);
                            }


                            b0pval = 2 * (1 - iagk0);
                            b1pval = 2 * (1 - iagk);
                            b2pval = 2 * (1 - iagk2);
                            b3pval = 2 * (1 - iagk3);
                            b4pval = 2 * (1 - iagk4);
                            b5pval = 2 * (1 - iagk5);
                            b6pval = 2 * (1 - iagk6);
                            b7pval = 2 * (1 - iagk7);
                            b8pval = 2 * (1 - iagk8);
                            b9pval = 2 * (1 - iagk9);
                            b10pval = 2 * (1 - iagk10);
                            b11pval = 2 * (1 - iagk11);
                            b12pval = 2 * (1 - iagk12);

                            iagkmatrix = new double[,] { { b0pval }, { b1pval }, { b2pval }, { b3pval }, { b4pval }, { b5pval }, { b6pval }, { b7pval }, { b8pval }, { b9pval }, { b10pval }, { b11pval }, { b12pval } };
                            w0 = new WaldTest(z0);
                            w1 = new WaldTest(z1);
                            w2 = new WaldTest(z2);
                            w3 = new WaldTest(z3);
                            w4 = new WaldTest(z4);
                            w5 = new WaldTest(z5);
                            w6 = new WaldTest(z6);
                            w7 = new WaldTest(z7);
                            w8 = new WaldTest(z8);
                            w9 = new WaldTest(z9);
                            w10 = new WaldTest(z10);
                            w11 = new WaldTest(z11);
                            w12 = new WaldTest(z12);
                            waldmatrix = new double[,] { { w0.PValue }, { w1.PValue }, { w2.PValue }, { w3.PValue }, { w4.PValue }, { w5.PValue }, { w6.PValue }, { w7.PValue }, { w8.PValue }, { w9.PValue }, { w10.PValue }, { w11.PValue }, { w12.PValue } };
                          
                        }
                       
                    }
                    catch (Exception ex)
                    {
                        Console.WriteLine("NaN's produced, " + dos.Columns[0] + " analysis did not converge");
                        goto SKIP2;
                    }

                    rwl.AcquireWriterLock(Timeout.Infinite);
                    if (print == "onlysnp" && names2.Count() > 1)
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            for (int b = 1; b < 2; b++)
                            {
                                if (int1 == null)
                                {
                                    jointint = -1;
                                }
                                if (chrome == null & file.Contains(".vcf"))
                                {
                                    chrome = splitFileContents2[0][0];
                                }
                                writer.WriteLine(names2[b] + " " + splitFileContents2[1][0] + " " + chrome + " " + splitFileContents2[2][0] + " " + splitFileContents2[3][0] + " " + splitFileContents2[4][0] + " " + Math.Round(NR[b, 0], 8) + " " + Math.Round(Math.Exp(NR[b, 0]), 8) + " " + Math.Round(Math.Exp(-NR[b, 0] * NR[names2.Count() - 1, 0]), 8) + " " + Math.Round(SE[b, b], 8) + " " + Math.Round(zscore[b, 0], 8) + " " + iagkmatrix[b, 0] + " " + waldmatrix[b, 0] + " " + jointint + " " + Math.Round(1 - allele, 8).ToString() + " " + Math.Round(allele, 8).ToString() + " " + Math.Round(infoscore, 8).ToString() + " " + Math.Round(NR[names2.Count() - 1, 0], 8));
                            }
                        }
                    }
                    else if (print == "onlyint" && names2.Count() > 1)
                    {

                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            for (int b = names2.Count() - 2; b < names2.Count() - 1; b++)
                            {
                                if (chrome == null & file.Contains(".vcf"))
                                {
                                    chrome = splitFileContents2[0][0];
                                }
                                writer.WriteLine(names2[b] + " " + splitFileContents2[1][0] + " " + chrome + " " + splitFileContents2[2][0] + " " + splitFileContents2[3][0] + " " + splitFileContents2[4][0] + " " + Math.Round(NR[b, 0], 8) + " " + Math.Round(Math.Exp(NR[b, 0]), 8) + " " + Math.Round(Math.Exp(-NR[b, 0] * NR[names2.Count() - 1, 0]), 8) + " " + Math.Round(SE[b, b], 8) + " " + Math.Round(zscore[b, 0], 8) + " " + iagkmatrix[b, 0] + " " + waldmatrix[b, 0] + " " + jointint + " " + Math.Round(1 - allele, 8).ToString() + " " + Math.Round(allele, 8).ToString() + " " + Math.Round(infoscore, 8).ToString() + " " + Math.Round(NR[names2.Count() - 1, 0], 8));
                            }
                        }
                    }
                    else
                    {
                        using (StreamWriter writer = new StreamWriter(SaveFile, true))
                        {
                            for (int b = 0; b < names2.Count(); b++)
                            {
                                if (int1 == null)
                                {
                                    jointint = -1;
                                }
                                if (chrome == null & file.Contains(".vcf"))
                                {
                                    chrome = splitFileContents2[0][0];
                                }
                                writer.WriteLine(names2[b] + " " + splitFileContents2[1][0] + " " + chrome + " " + splitFileContents2[2][0] + " " + splitFileContents2[3][0] + " " + splitFileContents2[4][0] + " " + Math.Round(NR[b, 0], 8) + " " + Math.Round(Math.Exp(NR[b, 0]), 8) + " " + Math.Round(Math.Exp(-NR[b, 0] * NR[names2.Count() - 1, 0]), 8) + " " + Math.Round(SE[b, b], 8) + " " + Math.Round(zscore[b, 0], 8) + " " + iagkmatrix[b, 0] + " " + waldmatrix[b, 0] + " " + jointint + " " + Math.Round(1 - allele, 8).ToString() + " " + Math.Round(allele, 8).ToString() + " " + Math.Round(infoscore, 8).ToString() + " " + Math.Round(NR[names2.Count() - 1, 0], 8));
                            }
                        }
                    }
                    rwl.ReleaseWriterLock();

                SKIP2:;
                    itercount = 0;
                    valid = 0;
                    H = null;
                    H1 = null;
                    V = null;
                    newArray = null;
                    newArray1 = null;
                    newArray2 = null;
                    newArray3 = null;
                    newArray4 = null;
                    newArray5 = null;
                    IV = null;
                    zscore = null;
                    iagkmatrix = null;
                    waldmatrix = null;
                    weiblike = null;
                    weiblike2 = null;
                    SE = null;
                    IH = null;
                    Z = null;
                    NR = null;
                    NRe = null;
                    ini = null;
                    ini1 = null;
                    ini2 = null;
                    Cov1 = null;
                    Cov2 = null;
                    Cov3 = null;
                    Cov4 = null;
                    Cov5 = null;
                    Cov6 = null;
                    Cov7 = null;
                    Cov8 = null;
                    Cov9 = null;
                    Cov10 = null;
                    SNP = null;
                    da = null;
                    db = null;
                    dc = null;
                    dd = null;
                    de = null;
                    df = null;
                    dg = null;
                    dh = null;
                    di = null;
                    dj = null;
                    dk = null;
                    dl = null;
                    d2a = null;
                    d2b = null;
                    d2c = null;
                    d2d = null;
                    d2e = null;
                    d2f = null;
                    d2g = null;
                    d2h = null;
                    d2i = null;
                    d2j = null;
                    d2k = null;
                    d2l = null;
                    dab = null;
                    dac = null;
                    dad = null;
                    dae = null;
                    daf = null;
                    dag = null;
                    dah = null;
                    dai = null;
                    daj = null;
                    dak = null;
                    dal = null;
                    dbc = null;
                    dbd = null;
                    dbe = null;
                    dbf = null;
                    dbg = null;
                    dbh = null;
                    dbi = null;
                    dbj = null;
                    dbk = null;
                    dbl = null;
                    dcd = null;
                    dce = null;
                    dcf = null;
                    dcg = null;
                    dch = null;
                    dci = null;
                    dcj = null;
                    dck = null;
                    dcl = null;
                    dde = null;
                    ddf = null;
                    ddg = null;
                    ddh = null;
                    ddi = null;
                    ddj = null;
                    ddk = null;
                    ddl = null;
                    def = null;
                    deg = null;
                    deh = null;
                    dei = null;
                    dej = null;
                    dek = null;
                    del = null;
                    dfg = null;
                    dfh = null;
                    dfi = null;
                    dfj = null;
                    dfk = null;
                    dfl = null;
                    dgh = null;
                    dgi = null;
                    dgj = null;
                    dgk = null;
                    dgl = null;
                    dhi = null;
                    dhj = null;
                    dhk = null;
                    dhl = null;
                    dij = null;
                    dik = null;
                    dil = null;
                    djk = null;
                    djl = null;
                    dkl = null;
                    dm = null;
                    d2m = null;
                    dma = null;
                    dmb = null;
                    dmc = null;
                    dmd = null;
                    dme = null;
                    dmf = null;
                    dmg = null;
                    dmh = null;
                    dmi = null;
                    dmj = null;
                    dmk = null;
                    dml = null;
                    w0 = null;
                    w1 = null;
                    w2 = null;
                    w3 = null;
                    w4 = null;
                    w5 = null;
                    w6 = null;
                    w7 = null;
                    w8 = null;
                    w9 = null;
                    w10 = null;
                    w11 = null;
                    w12 = null;
                    names2.Clear();
                }

            SKIP:;

                dos.Clear();
                dos.Columns.Clear();
                dos.Rows.Clear();
                independent.Clear();
                names.Clear();
                independent.Rows.Clear();
                independent.Columns.Clear();
                dt.Clear();
                dt.Columns.Clear();
                dt.Rows.Clear();
                dos = null;
                NN = null;
                PP = null;
                FF = null;
                FPsub = null;
                gsum2 = null;
                intColumn = null;
                independent = null;
                names = null;
                dt = null;
                pha = null;
                newnames = null;

            CRITFAIL:;

                table.Clear();
                table.Columns.Clear();
                table.Rows.Clear();
                table = null;
                splitFileContents2 = null;
                a = null;

                GC.Collect();

                counter++;
                randLineNum++;
            }
        FINISHED:;

            _semaSlim.Release();
            Console.WriteLine("Batch {0} analysis finished", threadId);

            if (_semaSlim.CurrentCount.Equals(int.Parse(cores)))
            {
                Console.WriteLine("Analysis complete!");
                Console.WriteLine("Congratulations " + Environment.UserName + ". You've earned the GWAS hero badge (level 1).");
            }

        }
        static void ShowHelp(OptionSet p)
        {
            p.WriteOptionDescriptions(Console.Out);
        }
    }
}

