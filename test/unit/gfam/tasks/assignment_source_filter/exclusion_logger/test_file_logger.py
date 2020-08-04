import unittest
from unittest.mock import mock_open, patch

from gfam.tasks.assignment_source_filter.exclusion_logger.file_logger import \
    FileLogger


class TestFileLogger(unittest.TestCase):

    @patch("builtins.open", new_callable=mock_open, read_data="data")
    def test_file_logger_should_open_file_and_write(self,
                                                    mock_file):
        # arrange
        file_name = "beautiful_logger_file"
        name = "protein_id"
        reason = "reason"
        expected_message = f"{name}: {reason}\n"

        # act
        _sut = FileLogger(file_name)
        _sut.log_exclusion(name, reason)

        # assert
        mock_file.assert_called_with(file_name, "a+")
        mock_file.return_value.__enter__().write.assert_called_once_with(
            expected_message)
