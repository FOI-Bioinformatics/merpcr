"""
Command-line interface for merPCR.
"""

import argparse
import logging
import sys
from typing import List

from .core.engine import (DEFAULT_IUPAC_MODE, DEFAULT_MARGIN,
                          DEFAULT_MISMATCHES, DEFAULT_PCR_SIZE,
                          DEFAULT_THREADS, DEFAULT_THREE_PRIME_MATCH,
                          DEFAULT_WORDSIZE, MerPCR)

# Default values for me-PCR compatibility
DEFAULT_MAX_STS_LINE_LENGTH = 1022


def convert_mepcr_arguments(args: List[str]) -> List[str]:
    """Convert me-PCR style arguments (M=50) to argparse style (-M 50)."""
    converted_args = []
    i = 0
    while i < len(args):
        arg = args[i]
        # Check for me-PCR style arguments (X=value)
        if len(arg) >= 3 and arg[1] == '=' and arg[0] in 'MNWXTQZISOP':
            param = arg[0]
            value = arg[2:]
            
            # Convert to argparse format
            if param == 'M':
                converted_args.extend(['-M', value])
            elif param == 'N':
                converted_args.extend(['-N', value])
            elif param == 'W':
                converted_args.extend(['-W', value])
            elif param == 'X':
                converted_args.extend(['-X', value])
            elif param == 'T':
                converted_args.extend(['-T', value])
            elif param == 'Q':
                converted_args.extend(['-Q', value])
            elif param == 'Z':
                converted_args.extend(['-Z', value])
            elif param == 'I':
                converted_args.extend(['-I', value])
            elif param == 'S':
                converted_args.extend(['-S', value])
            elif param == 'O':
                converted_args.extend(['-O', value])
            elif param == 'P':
                # P parameter is Mac-specific priority, ignore for now
                pass
        elif arg == '-help':
            # Convert me-PCR help to standard help
            converted_args.append('--help')
        else:
            # Keep all other arguments as is
            converted_args.append(arg)
        i += 1
    
    return converted_args


def setup_logging(quiet: int, debug: bool) -> None:
    """Set up logging based on arguments."""
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    logger = logging.getLogger("merpcr")

    if debug:
        logger.setLevel(logging.DEBUG)
    elif quiet == 0:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)


def margin_type(value):
    """Validate margin parameter."""
    ivalue = int(value)
    if ivalue < 0 or ivalue > 10000:
        raise argparse.ArgumentTypeError(f"Margin must be between 0-10000, got {ivalue}")
    return ivalue


def mismatch_type(value):
    """Validate mismatch parameter."""
    ivalue = int(value)
    if ivalue < 0 or ivalue > 10:
        raise argparse.ArgumentTypeError(f"Mismatches must be between 0-10, got {ivalue}")
    return ivalue


def wordsize_type(value):
    """Validate wordsize parameter."""
    ivalue = int(value)
    if ivalue < 3 or ivalue > 16:
        raise argparse.ArgumentTypeError(f"Word size must be between 3-16, got {ivalue}")
    return ivalue


def threads_type(value):
    """Validate threads parameter."""
    ivalue = int(value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError(f"Threads must be > 0, got {ivalue}")
    return ivalue


def pcr_size_type(value):
    """Validate PCR size parameter."""
    ivalue = int(value)
    if ivalue < 1 or ivalue > 10000:
        raise argparse.ArgumentTypeError(f"PCR size must be between 1-10000, got {ivalue}")
    return ivalue


def sts_line_length_type(value):
    """Validate STS line length parameter."""
    ivalue = int(value)
    if ivalue < 1:
        raise argparse.ArgumentTypeError(f"STS line length must be > 0, got {ivalue}")
    return ivalue


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser for merPCR."""
    parser = argparse.ArgumentParser(
        description="merPCR - Modern Electronic Rapid PCR",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("sts_file", type=str, help="STS file (tab-delimited)")
    parser.add_argument("fasta_file", type=str, help="FASTA sequence file")

    parser.add_argument(
        "-M",
        "--margin",
        type=margin_type,
        default=DEFAULT_MARGIN,
        help=f"Margin (default: {DEFAULT_MARGIN})",
    )

    parser.add_argument(
        "-N",
        "--mismatches",
        type=mismatch_type,
        default=DEFAULT_MISMATCHES,
        help=f"Number of mismatches allowed (default: {DEFAULT_MISMATCHES})",
    )

    parser.add_argument(
        "-W",
        "--wordsize",
        type=wordsize_type,
        default=DEFAULT_WORDSIZE,
        help=f"Word size (default: {DEFAULT_WORDSIZE})",
    )

    parser.add_argument(
        "-T",
        "--threads",
        type=threads_type,
        default=DEFAULT_THREADS,
        help=f"Number of threads (default: {DEFAULT_THREADS})",
    )

    parser.add_argument(
        "-X",
        "--three-prime-match",
        type=int,
        default=DEFAULT_THREE_PRIME_MATCH,
        help=f"Number of 3'-ward bases in which to disallow mismatches (default: {DEFAULT_THREE_PRIME_MATCH})",
    )

    parser.add_argument(
        "-O", "--output", type=str, default=None, help="Output file name (default: stdout)"
    )

    parser.add_argument(
        "-Q", "--quiet", type=int, choices=[0, 1], default=1, help="Quiet flag (0=verbose, 1=quiet)"
    )

    parser.add_argument(
        "-Z",
        "--default-pcr-size",
        type=pcr_size_type,
        default=DEFAULT_PCR_SIZE,
        help=f"Default PCR size (default: {DEFAULT_PCR_SIZE})",
    )

    parser.add_argument(
        "-I",
        "--iupac",
        type=int,
        choices=[0, 1],
        default=DEFAULT_IUPAC_MODE,
        help="IUPAC flag (0=don't honor IUPAC ambiguity symbols, 1=honor IUPAC symbols)",
    )

    parser.add_argument(
        "-S",
        "--max-sts-line-length",
        type=sts_line_length_type,
        default=DEFAULT_MAX_STS_LINE_LENGTH,
        help=f"Max. line length for the STS file (default: {DEFAULT_MAX_STS_LINE_LENGTH})",
    )

    parser.add_argument("-v", "--version", action="version", version="merPCR version 1.0.0")

    parser.add_argument("--debug", action="store_true", help="Enable debug logging")

    return parser


def main() -> int:
    """Main function to run the merPCR program."""
    # Convert me-PCR style arguments to argparse format
    converted_argv = convert_mepcr_arguments(sys.argv[1:])
    
    parser = create_parser()
    args = parser.parse_args(converted_argv)

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
            threads=args.threads,
            max_sts_line_length=args.max_sts_line_length,
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
