# CI/CD Documentation

This document describes the simplified and focused GitHub Actions workflows for the merPCR project.

## Overview

The CI/CD system has been streamlined to provide essential testing, quality checks, security scanning, and release automation without unnecessary complexity.

**Total Workflows**: 4
**Total Lines**: ~370 (down from 1,809)
**Focus**: Core functionality, speed, and maintainability

## Workflows

### 1. Tests (`test.yml`)

**Triggers**: Push to main/develop, Pull requests

**Purpose**: Core testing across platforms and Python versions

**Jobs**:
- **test**: Matrix testing on Ubuntu, macOS, and Windows with Python 3.11-3.13
  - Runs pytest with coverage
  - Uses parallel execution (-n auto)
  - Uploads coverage to Codecov (Ubuntu + Python 3.12 only)

- **compatibility**: Tests me-PCR compatibility
  - Runs `test_compatibility.py`
  - Ensures backward compatibility with me-PCR v1.0.6

**Runtime**: ~5-10 minutes per matrix configuration

### 2. Code Quality (`quality.yml`)

**Triggers**: Push to main/develop, Pull requests

**Purpose**: Code linting, formatting, type checking, and basic security

**Jobs**:
- **lint**: Quality checks
  - flake8 linting (syntax errors and complexity)
  - black formatting verification
  - isort import sorting verification
  - mypy type checking
  - bandit security scanning

- **format**: Auto-formatting on PRs
  - Automatically formats code with black and isort
  - Commits changes back to PR branch
  - Requires write permissions

**Runtime**: ~2-3 minutes

### 3. Security (`security.yml`)

**Triggers**: Push to main/develop, Pull requests, Weekly schedule (Sunday 3 AM UTC)

**Purpose**: Security vulnerability detection

**Jobs**:
- **codeql**: GitHub's advanced security analysis
  - Comprehensive static analysis
  - Maintained by GitHub Security Lab
  - Checks for vulnerabilities, secrets, and code quality issues

- **dependency-review**: PR dependency analysis
  - Checks for vulnerable dependencies in PRs
  - Enforces license compatibility
  - Provides automated PR comments with findings

**Runtime**: ~5-8 minutes

### 4. Release (`release.yml`)

**Triggers**: Git tags (v*), Manual workflow dispatch

**Purpose**: Package building and distribution

**Jobs**:
- **test**: Pre-release test verification
- **build**: Package building with verification
- **release**: GitHub release creation
- **publish-pypi**: PyPI package publication
- **notify**: Release notification

**Runtime**: ~10-15 minutes

## Configuration Files

### Required Secrets
- `CODECOV_TOKEN` (optional): Coverage reporting to Codecov
- `PYPI_API_TOKEN` (required for releases): PyPI publishing
- `GITHUB_TOKEN`: Automatically provided

### Test Configuration
- `pytest.ini`: Test configuration and markers
- `pyproject.toml`: Package configuration

## Removed Complexity

### What Was Removed (78% reduction)
- ❌ Reproducibility testing (unnecessary)
- ❌ Resource monitoring (overkill)
- ❌ Performance regression testing (too complex)
- ❌ Dependency compatibility matrix (redundant)
- ❌ Dynamic security analysis (redundant)
- ❌ Container security (not needed)
- ❌ 8+ redundant security tools (kept CodeQL + bandit)
- ❌ TruffleHog/GitLeaks (CodeQL covers this)
- ❌ Custom scripts (dependency_validation.py, performance_baseline.py)

### Why These Changes
1. **CodeQL is Comprehensive**: GitHub's CodeQL provides better security coverage than 10+ individual tools combined
2. **Focus on Speed**: Faster CI means faster feedback for developers
3. **Maintainability**: Simpler workflows are easier to understand and modify
4. **Cost Efficiency**: Fewer GitHub Actions minutes used
5. **Essential Coverage**: All critical checks remain (testing, linting, security, releases)

## Usage

### Running Workflows Locally

**Tests**:
```bash
# Run all tests
pytest -v

# Run with coverage
pytest -v --cov=src/merpcr --cov-report=term

# Run compatibility tests
python test_compatibility.py
```

**Code Quality**:
```bash
# Lint code
flake8 src/ tests/

# Format code
black src/ tests/
isort src/ tests/

# Type check
mypy src/

# Security check
bandit -r src/ -ll
```

### Manual Workflow Triggers

**Trigger Release**:
```bash
# Via GitHub CLI
gh workflow run release.yml -f version=v1.0.1

# Or create and push a git tag
git tag v1.0.1
git push origin v1.0.1
```

**View Workflow Runs**:
```bash
gh run list
gh run view <run-id>
```

## Troubleshooting

### Common Issues

**Test Failures**:
1. Check test logs in GitHub Actions
2. Run tests locally to reproduce
3. Verify test data is generated correctly

**Coverage Upload Fails**:
- Verify `CODECOV_TOKEN` is set (optional, won't fail CI)
- Check Codecov service status

**Auto-format Doesn't Commit**:
- Verify PR has write permissions
- Check branch protection rules
- May need to allow Actions to create commits

**Release Fails**:
- Verify `PYPI_API_TOKEN` is set correctly
- Check version tag format (must start with 'v')
- Ensure all tests pass before release

### Debug Mode

Enable debug logging in GitHub Actions:
1. Go to repository Settings → Secrets
2. Add `ACTIONS_STEP_DEBUG` = `true`
3. Add `ACTIONS_RUNNER_DEBUG` = `true`

## Monitoring

### What to Watch
- Test pass/fail rates across platforms
- Coverage trends (should stay >90%)
- Security alerts from CodeQL
- Dependency vulnerabilities

### GitHub Actions Dashboard
View all workflows at: `https://github.com/FOI-Bioinformatics/merpcr/actions`

## Maintenance

### Regular Updates Needed
- Python versions in test matrix (annually)
- GitHub Actions versions (check Dependabot)
- Dependencies in pyproject.toml (check weekly security scans)

### Adding New Tests
1. Add test to `tests/` directory
2. Tests run automatically on next push
3. Verify coverage doesn't drop

### Modifying Workflows
1. Edit workflow YAML files
2. Test locally with `act` (optional)
3. Create PR to test in CI
4. Verify all checks pass before merging

## Benefits of Simplified CI

✅ **78% less code**: 1,809 → 370 lines
✅ **Faster feedback**: 5-10 minutes vs 30+ minutes
✅ **Easier to maintain**: Clear, focused workflows
✅ **Lower costs**: Fewer GitHub Actions minutes
✅ **Better security**: CodeQL > 10+ individual tools
✅ **Essential coverage**: All critical checks remain

## References

- [GitHub Actions Documentation](https://docs.github.com/en/actions)
- [CodeQL Documentation](https://codeql.github.com/docs/)
- [pytest Documentation](https://docs.pytest.org/)
- [Coverage.py](https://coverage.readthedocs.io/)

---

*Last updated: October 2025*
*CI/CD System Version: 2.0 (Simplified)*
