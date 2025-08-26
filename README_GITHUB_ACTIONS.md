# GitHub Actions CI/CD Setup

This document explains the GitHub Actions workflows configured for the merPCR project.

## Workflows Overview

### 1. Continuous Integration (`ci.yml`)
**Triggers**: Push to main/develop branches, Pull requests
**Features**:
- Matrix testing across Python 3.9-3.12 on Ubuntu, macOS, Windows
- Automatic test data creation for missing files
- Test categorization: unit, integration, CLI, performance
- Coverage reporting with Codecov integration
- Performance benchmarking (on push/schedule only)

### 2. Code Quality (`code-quality.yml`)
**Triggers**: Push to main/develop branches, Pull requests
**Features**:
- Code linting with flake8
- Code formatting verification with black and isort
- Type checking with mypy (non-blocking)
- Security scanning with bandit
- Dependency vulnerability checking with safety
- Documentation validation
- Auto-formatting on pull requests (requires write permissions)

### 3. Security (`security.yml`)
**Triggers**: Push to main/develop, Pull requests, Weekly schedule
**Features**:
- Comprehensive security scanning with bandit, safety, pip-audit
- GitHub CodeQL analysis
- Dependency review for pull requests
- Secret detection with TruffleHog
- License compatibility checking

### 4. Release (`release.yml`)
**Triggers**: Git tags starting with 'v', Manual workflow dispatch
**Features**:
- Full test suite validation before release
- Package building and validation
- Automated GitHub releases
- PyPI publishing (requires PYPI_API_TOKEN secret)

## Configuration Files

### Test Configuration
- `pytest.ini`: pytest configuration with test markers
- `.coveragerc`: Coverage reporting settings
- `tox.ini`: Multi-environment testing configuration
- `Makefile`: Development convenience commands

### Issue Templates
- `.github/ISSUE_TEMPLATE/bug_report.yml`: Structured bug reports
- `.github/ISSUE_TEMPLATE/feature_request.yml`: Feature request template

## Fixed Issues

### Recent Fixes Applied:
1. **Deprecated Actions**: Updated to latest versions
   - `actions/upload-artifact@v3` → `v4`
   - `actions/download-artifact@v3` → `v4`
   - `actions/setup-python@v4` → `v5`
   - `github/codeql-action@v2` → `v3`

2. **Dependency Review**: Fixed configuration conflict
   - Removed `deny-licenses` to avoid conflict with `allow-licenses`
   - Updated to `actions/dependency-review-action@v4`

3. **Documentation Check**: Made non-blocking
   - Changed from error to warning for source newer than docs
   - Added error handling for missing files

4. **Auto-formatting Permissions**: Enhanced GitHub Actions permissions
   - Added explicit `contents: write` and `pull-requests: write` permissions
   - Updated checkout with proper fetch depth
   - Improved error handling for push failures

## Secrets Required

### Repository Secrets:
- `CODECOV_TOKEN` (optional): For coverage reporting
- `PYPI_API_TOKEN` (required for releases): PyPI publishing

### Automatic Secrets:
- `GITHUB_TOKEN`: Automatically provided by GitHub Actions

## Test Categories

The CI system recognizes these pytest markers:
- `unit`: Unit tests for individual components
- `integration`: End-to-end tests with real data
- `cli`: Command-line interface tests  
- `performance`: Performance and scalability benchmarks
- `slow`: Long-running tests (can be skipped)

## Performance Optimizations

### CI Environment Adaptations:
- Performance tests use smaller datasets in CI
- Automatic output redirection to prevent console spam
- Environment variable detection (`CI`, `GITHUB_ACTIONS`)
- Timeout protections for long-running tests

## Usage Examples

### Running Tests Locally:
```bash
# All tests
make test

# Specific categories
make test-unit
make test-integration
make test-performance

# With coverage
make coverage
```

### Manual Workflow Triggers:
```bash
# Trigger release workflow
gh workflow run release.yml -f version=v1.0.1

# View workflow runs
gh run list
```

## Troubleshooting

### Common Issues:
1. **Test failures**: Check test data creation in CI logs
2. **Permission errors**: Ensure repository has required secrets
3. **Coverage upload failures**: Check CODECOV_TOKEN configuration
4. **Release failures**: Verify PYPI_API_TOKEN and tag format

### Debug Steps:
1. Check workflow logs in GitHub Actions tab
2. Run tests locally with same environment
3. Verify all secrets are configured
4. Check branch protection rules

## Maintenance

### Regular Updates Needed:
- Python versions in matrix testing
- GitHub Actions versions
- Security scanning tool versions
- Dependencies in requirements

### Monitoring:
- Weekly security scans
- Performance trend analysis
- Test coverage changes
- Dependency updates

This CI/CD setup ensures code quality, security, and reliability across all supported platforms and Python versions.