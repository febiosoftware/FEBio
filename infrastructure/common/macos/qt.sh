#!/bin/zsh
set -e

REPO="https://github.com/qt/qtbase.git"
BRANCH="v6.6.1"

pushd $SOURCE_PATH
rm -rf "${SOURCE_PATH}/${BRANCH}"

git clone --depth 1 --branch $BRANCH $REPO $BRANCH

pushd $BRANCH
cmake . -G Ninja -L -B cmbuild \
  -DCMAKE_OSX_ARCHITECTURES="x86_64" \
  -DCMAKE_OSX_DEPLOYMENT_TARGET=10.15 \
  -DCMAKE_INSTALL_PREFIX="${INSTALLATION_PATH}/Qt" \
  -DCMAKE_PREFIX_PATH="${INSTALLATION_PATH}/homebrew" \
  -DFEATURE_developer_build=ON \
  -DFEATURE_private_tests=OFF \
  -DCMAKE_BUILD_TYPE=Release \
  -DQT_BUILD_TESTS=OFF \
  -DQT_BUILD_TESTS_BY_DEFAULT=OFF \
  -DQT_BUILD_EXAMPLES_BY_DEFAULT=OFF \
  -DQT_INSTALL_EXAMPLE_SOURCE_BY_DEFAULT=OFF \

pushd cmbuild
cmake --build . --parallel
cmake --install .
popd
popd

rm -rf "${BRANCH}"

cat << EOF >> /Users/$SSH_USER/.zprofile
export PATH="${INSTALLATION_PATH}/Qt/bin:\$PATH"
EOF

