#!/usr/bin/env sh
set -eu

IMAGE="${IMAGE:-gwswex-dev}"
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
CONTAINER_NAME="${CONTAINER_NAME:-gwswex-dev}"

usage() {
  cat <<'EOF'
Usage: scripts/dev-container.sh <command>

Commands:
  build      Build the image from .env.d/Containerfile
  run        Start a persistent container named gwswex-dev (for VS Code remote)
  shell      Open an interactive shell in /home with project bind-mounted
  test       Install package in editable mode and run pytest
  stop       Stop the persistent container
  remove     Remove the persistent container
EOF
}

detect_engine() {
  if command -v container >/dev/null 2>&1; then
    echo "container"
    return
  fi
  if command -v docker >/dev/null 2>&1; then
    echo "docker"
    return
  fi
  echo "Error: neither 'container' nor 'docker' is available on PATH." >&2
  exit 1
}

if [ "$#" -lt 1 ]; then
  usage
  exit 1
fi

ENGINE="$(detect_engine)"
CMD="$1"

case "$CMD" in
  build)
    "$ENGINE" build -f .env.d/Containerfile -t "$IMAGE" .
    ;;
  run)
    echo "Starting persistent container '$CONTAINER_NAME' from image '$IMAGE'..."
    if "$ENGINE" ps -a --filter "name=$CONTAINER_NAME" --format '{{.Names}}' | grep -q "^$CONTAINER_NAME$"; then
      echo "Container already exists. Removing old container..."
      "$ENGINE" rm -f "$CONTAINER_NAME"
    fi
    "$ENGINE" run -d --name "$CONTAINER_NAME" \
      --mount "type=bind,source=$PROJECT_DIR,target=/home" \
      -w /home \
      "$IMAGE" \
      sh -c "while true; do sleep 1; done"
    echo "Container '$CONTAINER_NAME' is running. You can connect VS Code to it now."
    ;;
  shell)
    "$ENGINE" run --rm -it \
      --mount "type=bind,source=$PROJECT_DIR,target=/home" \
      -w /home \
      "$IMAGE"
    ;;
  test)
    "$ENGINE" run --rm \
      --mount "type=bind,source=$PROJECT_DIR,target=/home" \
      -w /home \
      "$IMAGE" \
      sh -lc "pip install -e . --no-build-isolation && pytest tests/"
    ;;
  stop)
    echo "Stopping container '$CONTAINER_NAME'..."
    "$ENGINE" stop "$CONTAINER_NAME" || echo "Container not running."
    ;;
  remove)
    echo "Removing container '$CONTAINER_NAME'..."
    "$ENGINE" rm -f "$CONTAINER_NAME" || echo "Container does not exist."
    ;;
  *)
    usage
    exit 1
    ;;
esac
