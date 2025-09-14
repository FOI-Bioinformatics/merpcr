/**
 * @name Bioinformatics Security Analysis
 * @description Security patterns specific to bioinformatics and genomic data processing
 * @kind problem
 * @problem.severity warning
 * @precision medium
 * @id python/bioinformatics-security
 * @tags security
 *       bioinformatics
 *       genomics
 *       file-handling
 */

import python
import semmle.python.security.dataflow.UnsafeShellQuery

/**
 * Detect potential file path traversal vulnerabilities in genomic file handling
 */
class UnsafeGenomicFileHandling extends TaintTracking::Configuration {
  UnsafeGenomicFileHandling() { this = "UnsafeGenomicFileHandling" }

  override predicate isSource(DataFlow::Node node) {
    exists(Call call |
      call.getFunc().(Name).getId() in ["input", "raw_input"] and
      node.asExpr() = call
    )
  }

  override predicate isSink(DataFlow::Node node) {
    exists(Call call |
      call.getFunc().(Attribute).getName() = "open" and
      node.asExpr() = call.getArg(0)
    )
    or
    exists(Call call |
      call.getFunc().(Name).getId() in ["load_sts_file", "load_fasta_file"] and
      node.asExpr() = call.getArg(0)
    )
  }
}

/**
 * Detect hardcoded genomic file paths that might contain sensitive data
 */
predicate hardcodedGenomicPath(StrConst path) {
  path.getText().regexpMatch(".*(home|users)/[^/]+/(.*\\.(fa|fasta|fastq|fq|sts|sam|bam|vcf|gff|gtf)).*")
}

/**
 * Detect potential buffer overflow in sequence processing
 */
predicate unsafeSequenceProcessing(Call call) {
  call.getFunc().(Attribute).getName() in ["read", "readline", "readlines"] and
  not exists(Call sizeLimit |
    sizeLimit.getFunc().(Attribute).getName() = "read" and
    sizeLimit.getNumArgs() > 0
  )
}

/**
 * Detect missing input validation for genomic data formats
 */
predicate missingGenomicValidation(FunctionDef func) {
  func.getName().regexpMatch(".*(load|parse|read).*(sts|fasta|fastq).*") and
  not exists(If validation |
    validation.getParent*() = func and
    validation.getTest().(Compare).getLeft().(Call).getFunc().(Attribute).getName() in 
      ["exists", "isfile", "isdir", "endswith", "startswith"]
  )
}

/**
 * Results
 */
from Expr problem, string message
where
  (
    exists(UnsafeGenomicFileHandling config, DataFlow::PathNode source, DataFlow::PathNode sink |
      config.hasFlowPath(source, sink) and
      problem = sink.getNode().asExpr() and
      message = "Potential path traversal vulnerability in genomic file handling"
    )
  ) or
  (
    hardcodedGenomicPath(problem) and
    message = "Hardcoded genomic file path may expose sensitive data location"
  ) or
  (
    unsafeSequenceProcessing(problem) and
    message = "Potentially unsafe genomic sequence processing without size limits"
  ) or
  (
    exists(FunctionDef func |
      missingGenomicValidation(func) and
      problem = func and
      message = "Genomic data loading function lacks input validation"
    )
  )
select problem, message