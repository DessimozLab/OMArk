'''
            OMArk - Quality assesment of coding-gene repertoire annotation
            (C) 2022 Yannis Nevers <yannis.nevers@unil.ch>
            This file is part of OMArk.
            OMArk is free software: you can redistribute it and/or modify
            it under the terms of the GNU Lesser General Public License as published by
            the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.
            OMArk is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
            GNU Lesser General Public License for more details.
            You should have received a copy of the GNU Lesser General Public License
            along with OMArk. If not, see <http://www.gnu.org/licenses/>.
        '''
import logging
import sys

logging.root.handlers = []
h1 = logging.StreamHandler(sys.stdout)
h1.addFilter(lambda record: record.levelno <= logging.INFO)
h2 = logging.StreamHandler(sys.stderr)
h2.addFilter(lambda record: record.levelno > logging.INFO)
logging.basicConfig(format="%(levelname)s: %(message)s", handlers  = [h1,h2])
LOG = logging.getLogger(__name__)


def set_log_level(x):
    LOG.setLevel(getattr(logging, x.upper()))
