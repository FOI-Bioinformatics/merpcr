# merPCR Documentation

Welcome to the merPCR documentation. This directory contains comprehensive guides and reference materials for using merPCR effectively.

## Documentation Structure

### 📖 [User Guide](USER_GUIDE.md)
Comprehensive guide for using merPCR, including:
- Getting started and installation
- Understanding STS markers
- Data preparation
- Parameter tuning
- Performance optimization
- Troubleshooting

### 📚 [API Reference](API.md) 
Detailed technical documentation including:
- Class and method specifications
- Parameter descriptions
- Return values and exceptions
- Internal method documentation
- Compatibility notes

### 💡 [Examples and Tutorials](EXAMPLES.md)
Practical examples and step-by-step tutorials:
- Basic usage examples
- Real-world case studies
- Python API examples
- Batch processing workflows
- Advanced analysis techniques

## Quick Navigation

### For New Users
1. Start with the [User Guide](USER_GUIDE.md#getting-started)
2. Follow [Tutorial 1](EXAMPLES.md#tutorial-1-first-time-setup) for hands-on experience
3. Review [basic examples](EXAMPLES.md#basic-examples) for common use cases

### For Developers
1. Review the [API Reference](API.md) for technical details
2. Check [Python API examples](EXAMPLES.md#python-api-examples) for programmatic usage
3. See [advanced examples](EXAMPLES.md#advanced-examples) for complex scenarios

### For Specific Use Cases
- **Genome analysis**: [Tutorial 2](EXAMPLES.md#tutorial-2-human-genome-analysis)
- **Cross-species comparison**: [Tutorial 3](EXAMPLES.md#tutorial-3-cross-species-comparison)
- **High-throughput processing**: [Tutorial 4](EXAMPLES.md#tutorial-4-high-throughput-analysis)
- **Performance optimization**: [Performance Guide](USER_GUIDE.md#performance-optimization)

## Common Tasks

### Installing merPCR
See [Installation Guide](USER_GUIDE.md#getting-started)

### Running Your First Analysis
```bash
# Basic command
merpcr markers.sts genome.fa

# With common options
merpcr -M 50 -N 1 -T 4 -O results.txt markers.sts genome.fa
```

### Understanding Parameters
See [Parameter Tuning](USER_GUIDE.md#parameter-tuning) and [API Reference](API.md#constructor)

### Interpreting Results
See [Interpreting Results](USER_GUIDE.md#interpreting-results)

### Troubleshooting Issues
See [Troubleshooting](USER_GUIDE.md#troubleshooting)

## Additional Resources

### Testing
- Run `make test` for comprehensive testing
- See testing documentation in the main [README](../README.md#testing)

### Performance
- Use `make test-performance` for benchmarking
- See [Performance Optimization](USER_GUIDE.md#performance-optimization)

### Contributing
- Follow the development setup in the main [README](../README.md)
- Review the [API Reference](API.md) for implementation details

## Support

If you encounter issues or have questions:

1. Check the [Troubleshooting](USER_GUIDE.md#troubleshooting) section
2. Review relevant [examples](EXAMPLES.md) for similar use cases  
3. Consult the [API Reference](API.md) for technical details
4. Enable debug mode (`--debug`) for detailed error information

## Document Versions

These documents are maintained alongside the merPCR codebase and reflect the current version. For historical versions, check the git repository history.

Last updated: 2024 (Version 1.0.0)