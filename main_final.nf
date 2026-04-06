
// main.simple.v2.nf — Minimal STAR + RSEM orchestrator (DSL2), no top-level publishDir
nextflow.enable.dsl=2

/*
  Requirements:
    - Pass explicit references:
        --genomeDir  /path/to/STAR/index
        --rsem_ref   /path/to/RSEM/prefix
    - Provide input as either:
        --samplesheet path/to/sheet.csv
      OR
        --reads '/path/*_{R1,R2}.fastq.gz'   (use --pe false for single-end)
*/

// ---------------- Params (minimal, with sane defaults) ----------------
params.outdir        = params.outdir        ?: 'results'
params.strandedness  = params.strandedness  ?: 'unstranded'  // or reverse if dUPT
params.genomeDir     = params.genomeDir     ?: null     // REQUIRED
params.rsem_ref      = params.rsem_ref      ?: null     // REQUIRED (prefix string)
params.samplesheet   = params.samplesheet   ?: null
params.rscript       = params.rscript       ?: null
params.reads         = params.reads         ?: null
params.pe            = (params.pe in [false,'false',0,'0']) ? false : true
params.stranded      = params.stranded      ?: 'none'   // none|forward|reverse
params.read_group    = params.read_group    ?: null     // optional RG text (module-specific)
params.read_cmd      = params.read_cmd      ?: 'zcat'   // zcat for .gz
params.t1k_hla_fasta = params.t1k_hla_fasta ?: null   // REQUIRED for T1K (e.g. .../hlaidx_rna_seq.fa)


// ---------------- Modules ----------------
include { STAR_ALIGN ; RSEM_CALCULATE  } from './modules/star_rsem.nf'
include { samtools_index } from './modules/util.nf'
include { PREP_STAR_INDEX } from './modules/indexing.nf'

//include { RSEM_POSTREPORT_MIN } from './modules/make_report.nf'

//include { RMATS_single_sample } from './modules/rmats.nf'
include { FASTQC} from './modules/QC.nf'
include { T1K_TYPING } from './modules/t1k_hla.nf'
//include { RSCRIPT_PLOTS } from './modules/star_rsem_plot.nf'

// ---------------- Helpers ----------------
def inferPairGroup(File f) {
  def n = f.getName()
  n = n
    .replaceAll(/(?i)_R[12]\b/, '')
    .replaceAll(/-[12]\b/, '')
    .replaceAll(/(?:_1|_2)\b/, '')
    .replaceAll(/(?i)\.fastq(\.gz)?$/, '')
  return n
}

// ---------------- Validate references early ----------------
if( !params.genomeDir ) error "Missing --genomeDir (STAR index directory)"
if( !file(params.genomeDir).exists() ) error "STAR index not found: ${params.genomeDir}"
if( !params.rsem_ref ) error "Missing --rsem_ref (RSEM reference prefix, e.g. /refs/mouse/.../gencode_vM33)"
if( !params.t1k_hla_fasta ) error "Missing --t1k_hla_fasta (e.g. /path/to/hlaidx/hlaidx_rna_seq.fa)"
if( !file(params.t1k_hla_fasta).exists() ) error "T1K HLA fasta not found: ${params.t1k_hla_fasta}"
if( !params.pe ) error "T1K_TYPING is configured for paired-end only. Run with --pe true (default) and provide R1/R2."




// ---------------- Build input samples channel ----------------
def samples_ch

if( params.samplesheet ) {
  log.info "Loading samplesheet: ${params.samplesheet}"
  samples_ch = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row ->
      // Accept both header styles: (sample,fastq1,fastq2,read_group|strandedness)
      def sid = (row.sample ?: row.sample_id) as String
      def f1  = (row.fastq1 ?: row.fastq_1)?.toString()?.trim()
      def f2  = (row.fastq2 ?: row.fastq_2)?.toString()?.trim()
      def r1  = file(f1)
      def r2  = f2 ? file(f2) : null

      if( !r1.exists() ) error "FASTQ not found for sample ${sid}: ${r1}"
      if( params.pe && !r2 ) error "Paired-end requested but fastq2 missing for sample ${sid}"

      // read group: prefer column, then param, else synthesize
      def rg_raw = (row.read_group ?: params.read_group)
      def rg     = rg_raw ? rg_raw.toString() : "ID:${sid};SM:${sid}"

      tuple(sid, r1, r2, rg)
    }
} else if( params.reads ) {
  log.info "Globbing reads: ${params.reads}"
  if( params.pe ) {
    samples_ch = Channel
      .fromPath(params.reads)
      .groupBy { inferPairGroup(it) }
      .map { group, files ->
        def isR1 = { s -> s ==~ /(?i).*(?:_R1\b|\bR1\b|_1\b|-1\b).*/ }
        def isR2 = { s -> s ==~ /(?i).*(?:_R2\b|\bR2\b|_2\b|-2\b).*/ }
        def r1 = files.find { f -> isR1(f.name) }
        def r2 = files.find { f -> isR2(f.name) }
        if( !r1 || !r2 ) error "Could not auto-pair reads for group ${group}. Found: ${files*.name}"
        def rg = params.read_group ?: "ID:${group};SM:${group}"
        tuple(group as String, r1, r2, rg)
      }
  } else {
    samples_ch = Channel
      .fromPath(params.reads)
      .map { f ->
        def rg = params.read_group ?: "ID:${f.baseName};SM:${f.baseName}"
        tuple(f.baseName, f, null, rg)
      }
  }
} else {
  error "Provide either --samplesheet or --reads"
}

// ---------------- Constant reference channels ----------------
def star_index_ch  = Channel.value( file(params.genomeDir) )
def rsem_prefix_ch = Channel.value( params.rsem_ref as String )

// ---------------- Normalize tuple channel to avoid Broadcast cast ----------------
samples_ch = samples_ch.map { it }


// ---------------- Workflow ----------------
workflow {
  FASTQC(samples_ch)
  samples_ch.view { "SAMPLES_CH item => ${it}" }
  // Broadcast samples to multiple consumers
  //samples_ch.into { samples_fastqc_ch; samples_star_ch; samples_t1k_ch }



  aligned    = STAR_ALIGN( samples_ch, star_index_ch )
  T1K_TYPING( samples_ch.map { sid, r1, r2, rg -> tuple(sid, r1, r2) } )

  samtools_index(aligned.bam)
  quantified = RSEM_CALCULATE( aligned.tx_bam.map{ sid, bam -> tuple(sid,bam) }, rsem_prefix_ch )


  //rscript_plots(samples_ch)
 
   // run the subworkflow
  //RSEM_POSTREPORT_MIN( quantified )

  // collect its named outputs
  //report_pdf_ch = RSEM_POSTREPORT_MIN.out.report_pdf
  //summaries_ch  = RSEM_POSTREPORT_MIN.out.summaries


  // optional: publish or use downstream
  //report_pdf_ch, summaries_ch = RSEM_POSTREPORT_MIN( quantified )

  // optional: peek
  //report_pdf_ch.view { "RSEM report -> ${it}" }
   
}
