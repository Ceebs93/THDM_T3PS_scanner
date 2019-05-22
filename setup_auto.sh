#!/usr/bin/env bash

echo -e "source env.sh"
source env.sh

echo -e "./packages/install.sh"
./packages/install.sh


echo -e "./setup_links.sh"
./setup_links.sh
