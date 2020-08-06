import unittest
from unittest.mock import patch, mock_open

from gfam.tasks.assignment_source_filter.exclusion_logger.empty_logger import \
    EmptyLogger
from gfam.tasks.assignment_source_filter.exclusion_logger.exclusion_logger\
    import ExclusionLogger
from gfam.tasks.assignment_source_filter.exclusion_logger.file_logger import \
    FileLogger


class TestExclusionLogger(unittest.TestCase):
    def test_logger_without_file_creates_empty_logger(self):
        # arrange
        # act
        _sut = ExclusionLogger()

        # assert
        self.assertIsInstance(_sut._exclusion_logger, EmptyLogger)

    @patch("builtins.open", new_callable=mock_open, read_data="data")
    def test_logger_with_file_creates_file_logger_and_opens_file(self,
                                                                 mock_file):
        # arrange
        file_name = "beautiful_logger_file"

        # act
        _sut = ExclusionLogger(file_name)

        # assert
        self.assertIsInstance(_sut._exclusion_logger, FileLogger)
        mock_file.assert_called_with(file_name, "a+")
