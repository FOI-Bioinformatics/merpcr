"""
Command-line interface for merPCR.
"""

import argparse
import logging
import sys

from .core.engine import MerPCR, DEFAULT_MARGIN, DEFAULT_WORDSIZE, DEFAULT_MISMATCHES, \
    DEFAULT_THREE_PRIME_MATCH, DEFAULT_PCR_SIZE, DEFAULT_IUPAC_MODE, DEFAULT_THREADS


def setup_logging(quiet: int, debug: bool) -> None:
    """Set up logging based on arguments."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    
    logger = logging.getLogger("merpcr")
    
    if debug:
        logger.setLevel(logging.DEBUG)
    elif quiet == 0:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser for merPCR."""
    parser = argparse.ArgumentParser(
        description="merPCR - Modern Electronic Rapid PCR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('sts_file', type=str, help="STS file (tab-delimited)")
    parser.add_argument('fasta_file', type=str, help="FASTA sequence file")

    parser.add_argument('-M', '--margin', type=int, default=DEFAULT_MARGIN,
                        help=f"Margin (default: {DEFAULT_MARGIN})")

    parser.add_argument('-N', '--mismatches', type=int, default=DEFAULT_MISMATCHES,
                        help=f"Number of mismatches allowed (default: {DEFAULT_MISMATCHES})")

    parser.add_argument('-W', '--wordsize', type=int, default=DEFAULT_WORDSIZE,
                        help=f"Word size (default: {DEFAULT_WORDSIZE})")

    parser.add_argument('-T', '--threads', type=int, default=DEFAULT_THREADS,
                        help=f"Number of threads (default: {DEFAULT_THREADS})")

    parser.add_argument('-X', '--three-prime-match', type=int, default=DEFAULT_THREE_PRIME_MATCH,
                        help=f"Number of 3'-ward bases in which to disallow mismatches (default: {DEFAULT_THREE_PRIME_MATCH})")

    parser.add_argument('-O', '--output', type=str, default=None,
                        help="Output file name (default: stdout)")

    parser.add_argument('-Q', '--quiet', type=int, choices=[0, 1], default=1,
                        help="Quiet flag (0=verbose, 1=quiet)")

    parser.add_argument('-Z', '--default-pcr-size', type=int, default=DEFAULT_PCR_SIZE,
                        help=f"Default PCR size (default: {DEFAULT_PCR_SIZE})")

    parser.add_argument('-I', '--iupac', type=int, choices=[0, 1], default=DEFAULT_IUPAC_MODE,
                        help="IUPAC flag (0=don't honor IUPAC ambiguity symbols, 1=honor IUPAC symbols)")

    parser.add_argument('-v', '--version', action='version',
                        version="merPCR version 1.0.0")

    parser.add_argument('--debug', action='store_true',
                        help="Enable debug logging")

    return parser


def main() -> int:
    """Main function to run the merPCR program."""
    parser = create_parser()
    args = parser.parse_args()

    # Set up logging
    setup_logging(args.quiet, args.debug)
    
    logger = logging.getLogger("merpcr")

    try:
        # Initialize the merPCR instance with command line arguments
        mer_pcr = MerPCR(
            wordsize=args.wordsize,
            margin=args.margin,
            mismatches=args.mismatches,
            three_prime_match=args.three_prime_match,
            iupac_mode=args.iupac,
            default_pcr_size=args.default_pcr_size,
            threads=args.threads
        )

        # Load STS file
        if not mer_pcr.load_sts_file(args.sts_file):
            logger.error(f"Failed to load STS file: {args.sts_file}")
            return 1

        # Load FASTA file
        fasta_records = mer_pcr.load_fasta_file(args.fasta_file)
        if not fasta_records:
            logger.error(f"Failed to load FASTA file: {args.fasta_file}")
            return 1

        # Run the search
        hit_count = mer_pcr.search(fasta_records, args.output)

        logger.info(f"Search complete: {hit_count} hits found")
        return 0

    except Exception as e:
        logger.error(f"Error: {str(e)}")
        if args.debug:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())