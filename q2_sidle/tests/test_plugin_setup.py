from unittest import TestCase, main

from q2_sidle.plugin_setup import plugin

class PluginSetupTest(TestCase):
	def test_plugin_setup(self):
		self.assertEqual(plugin.name, 'sidle')

if __name__ == '__main__':
    main()