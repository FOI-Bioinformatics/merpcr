# Continuous Integration and Deployment Infrastructure

This document provides technical specifications for the GitHub Actions automation workflows implemented for the merPCR project.

## Workflow Architecture

### 1. Continuous Integration Pipeline (`ci.yml`)
**Activation Triggers**: Repository push events to main/develop branches, Pull request submissions
**Technical Capabilities**:
- Cross-platform validation matrix spanning Python 3.8-3.12 across Ubuntu, macOS, and Windows environments
- Automated test dataset generation for missing reference files
- Systematic test categorization: unit validation, integration testing, CLI verification, performance benchmarking
- Code coverage analysis with Codecov reporting integration
- Computational performance regression testing (activated on push/scheduled intervals)

### 2. Code Quality Assurance (`code-quality.yml`)
**Activation Triggers**: Repository modifications to main/develop branches, Pull request events
**Quality Control Features**:
- Static code analysis via flake8 linting framework
- Code standardization verification using black and isort formatters
- Type annotation validation through mypy analysis (advisory mode)
- Security vulnerability assessment via bandit scanning
- Dependency security evaluation using safety auditing tools
- Documentation consistency validation
- Automated code formatting on pull request submissions (requires repository write permissions)

### 3. Security Assessment Framework (`security.yml`)
**Activation Triggers**: Repository updates to main/develop branches, Pull request submissions, Weekly automated schedule
**Security Capabilities**:
- Multi-layered security analysis utilizing bandit, safety, and pip-audit tools
- GitHub CodeQL semantic code analysis
- Dependency vulnerability assessment for pull request modifications
- Credential and secret detection through TruffleHog scanning
- Software license compatibility verification

### 4. Release Management Pipeline (`release.yml`)
**Activation Triggers**: Git tag creation with 'v' prefix, Manual workflow dispatch
**Release Features**:
- Comprehensive test suite validation prior to package creation
- Distribution package construction and integrity validation
- Automated GitHub release generation with metadata
- PyPI distribution publishing (requires PYPI_API_TOKEN authentication)

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

## Test Classification Framework

The continuous integration system implements the following pytest marker taxonomy:
- `unit`: Isolated component validation tests
- `integration`: End-to-end workflow validation using authentic datasets
- `cli`: Command-line interface functionality verification
- `performance`: Computational efficiency and scalability benchmarking
- `slow`: Extended-duration tests (optionally excluded for rapid validation)

## Computational Optimizations

### CI Environment Adaptations:
- Performance validation utilizes reduced datasets for CI environment compatibility
- Automated output management to minimize console verbosity
- Environment variable detection for CI/CD context awareness (`CI`, `GITHUB_ACTIONS`)
- Execution timeout safeguards for extended-duration test procedures

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

### Common Diagnostic Issues:
1. **Test execution failures**: Examine test dataset generation procedures in CI execution logs
2. **Authentication failures**: Confirm repository secret configuration completeness
3. **Coverage reporting failures**: Validate CODECOV_TOKEN authentication configuration
4. **Release deployment failures**: Verify PYPI_API_TOKEN credentials and tag format compliance

### Systematic Debugging Protocol:
1. Examine detailed workflow execution logs via GitHub Actions interface
2. Replicate test execution in local environment matching CI specifications
3. Confirm complete authentication secret configuration
4. Validate branch protection rule compliance

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

This comprehensive CI/CD infrastructure ensures systematic code quality assurance, security validation, and operational reliability across all supported computational platforms and Python version specifications.