import unittest
from unittest.mock import MagicMock, call, patch

from gfam.tasks.assignment_source_filter.exclusion_logger.file_logger import \
    FileLogger


class TestFileLogger(unittest.TestCase):

    @patch("builtins.open", new_callable=MagicMock())
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
        _sut.close()

        # assert
        mock_file.assert_called_with(file_name, "a+")
        self.assertIn(call().write(expected_message), mock_file.mock_calls)
        self.assertIn(call().close(), mock_file.mock_calls)
