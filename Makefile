# Makefile for merPCR

.PHONY: all test clean install docs

# Default target
all: test

# Install the package
install:
	pip install -e .

# Run tests
test:
	python -m unittest test_merPCR.py

# Run a simple demonstration with the test data
demo:
	python merPCR.py test/test.sts test/test.fa

# Clean up temporary files
clean:
	rm -rf __pycache__
	rm -rf *.egg-info
	rm -rf build/
	rm -rf dist/
	rm -f *.pyc
	find . -name "__pycache__" -type d -exec rm -rf {} +
	find . -name "*.pyc" -delete

# Format code with Black (must be installed)
format:
	black merPCR.py test_merPCR.py example.py setup.py

# Run static type checking with mypy (must be installed)
typecheck:
	mypy merPCR.py

# Generate HTML documentation with pydoc
docs:
	mkdir -p docs
	pydoc -w merPCR
	mv merPCR.html docs/