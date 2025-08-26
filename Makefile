.PHONY: test coverage lint format clean install dev-install build upload help

help:
	@echo "Available targets:"
	@echo "  test          - Run all tests"
	@echo "  test-unit     - Run unit tests only"
	@echo "  test-integration - Run integration tests only"
	@echo "  test-performance - Run performance tests"
	@echo "  coverage      - Run tests with coverage report"
	@echo "  lint          - Run code linting"
	@echo "  format        - Format code with black and isort"
	@echo "  clean         - Clean build artifacts"
	@echo "  install       - Install package"
	@echo "  dev-install   - Install package in development mode"
	@echo "  build         - Build distribution packages"
	@echo "  upload        - Upload to PyPI (requires authentication)"

test:
	pytest

test-unit:
	pytest -m "unit or not (integration or performance or cli)"

test-integration:
	pytest -m integration

test-performance:
	SKIP_PERFORMANCE_TESTS= pytest tests/test_performance.py -v

coverage:
	coverage run -m pytest
	coverage report
	coverage html
	@echo "Coverage report generated in htmlcov/index.html"

lint:
	flake8 src/ tests/
	black --check src/ tests/
	isort --check-only src/ tests/

format:
	black src/ tests/
	isort src/ tests/

clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .tox/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf coverage.xml
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

install:
	pip install .

dev-install:
	pip install -e .

build: clean
	python -m build

upload: build
	python -m twine upload dist/*