from django.test import TestCase, Client
from django.urls import reverse


class ExploreViewTest(TestCase):
    def test_main_explore_view(self):
        # Create a client
        client = Client()

        # get request to the explore view
        url = reverse("explore:index")
        response = client.get(url)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "explore/main_explore.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "explore")
