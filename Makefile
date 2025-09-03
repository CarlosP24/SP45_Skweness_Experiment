CLUSTERS_CFG ?= bin/clusters.yaml
CLUSTER ?= esbirro
ARGUMENT ?=

SHORT := $(shell yq e '.${CLUSTER}.short' $(CLUSTERS_CFG))
HOST := $(shell yq e '.${CLUSTER}.host' $(CLUSTERS_CFG))
RPATH1 := $(shell yq e '.${CLUSTER}.path' $(CLUSTERS_CFG))
CUR_DIR := $(notdir $(CURDIR))
RPATH := $(RPATH1)/$(CUR_DIR)

LAUNCHER := $(shell yq e '.${CLUSTER}.launcher' $(CLUSTERS_CFG))

SRC_DIR := .
EXCLUDE := --exclude='.git'
RSYNC_OPTS := -ax --delete $(EXCLUDE)

# ----------------------------------------
# Phony Targets
# ----------------------------------------
.PHONY: all deploy run clean

# ----------------------------------------
# Default Target
# ----------------------------------------
all:
	@echo "Available targets:"
	@echo "  make deploy       # Sync local code to $(CLUSTER)"
	@echo "  make run     # Run tests on $(CLUSTER)"
	@echo "  make clean        # Clean local build artifacts"
	@echo ""
	@echo "Override CLUSTER on the command line:"
	@echo "  make CLUSTER=clusterB deploy"


# ----------------------------------------
# Deploy: sync code to remote cluster
# ----------------------------------------
deploy:
	@if [ -z "$(SHORT)" ]; then \
	  echo "Error: Unknown cluster '$(CLUSTER)' in $(CLUSTERS_CFG)"; \
	  exit 1; \
	fi
	@echo "Deploying to $(CLUSTER) @ $(SHORT):$(RPATH)"
	ssh $(SHORT) "mkdir -p $(RPATH)"
	rsync $(RSYNC_OPTS) $(SRC_DIR)/ $(SHORT):$(RPATH)

# ----------------------------------------
# Test-Run: execute test suite on remote
# ----------------------------------------
run: deploy
	@echo "Running on $(CLUSTER) @ $(SHORT)"
	ssh -A $(SHORT) "bash -l -c 'cd $(RPATH) && bash $(LAUNCHER) $(ARGUMENT)'"

# ----------------------------------------
# Clean: remove local build artifacts
# ----------------------------------------
clean:
	@echo "Cleaning local build artifacts..."
	@find . -type f \( -name '*.o' -o -name '*.pyc' \) -delete
	@rm -rf build/