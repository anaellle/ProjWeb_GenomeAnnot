from django.test import TestCase, Client
from django.urls import reverse

from main.models import CustomUser

class ViewTest(TestCase):

    def test_main_annotate_view_base(self):
        # Create a client
        self.client = Client()
        # get request to the annotate view
        url = reverse("main:annotate")
        response = self.client.get(url,follow=True) #follow to go to each redirection
        # Check 200 status code
        self.assertEqual(response.status_code, 200)
        # Check template
        self.assertTemplateUsed(response, "main/home.html")
        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "home")

    def test_main_annotate_view(self):
        # Créer un super utilisateur
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)
        # get request to the annotate view
        url = reverse("main:annotate")
        response = client.get(url, follow=True)
        # Check 200 status code
        self.assertEqual(response.status_code, 200)
        # Check template (décommentez cette ligne si vous souhaitez vérifier le modèle utilisé)
        self.assertTemplateUsed(response, "main/annotate/main_annotate.html")
        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "annotate")

    def test_main_blast_view(self):
        # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)
        # get request to the blast view
        url = reverse("main:blast")
        response = client.get(url,follow=True)
        # Check 200 status code
        self.assertEqual(response.status_code, 200)
        # Check template
        self.assertTemplateUsed(response, "main/blast/main_blast.html")
        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "blast")

    def test_main_explore_view(self):
        # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)
        # get request to the explore view
        url = reverse("main:explore")
        response = client.get(url,follow=True)
        # Check 200 status code
        self.assertEqual(response.status_code, 200)
        # Check template
        self.assertTemplateUsed(response, "main/explore/main_explore.html")
        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "explore")

    def test_main_validate_view(self):
        # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)
        # get request to the validate view
        url = reverse("main:validate")
        response = client.get(url,follow=True)
        # Check 200 status code
        self.assertEqual(response.status_code, 200)
        # Check template
        self.assertTemplateUsed(response, "main/validate/main_validate.html")
        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "validate")

    def test_home_view(self):
       # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)

        # get request to the validate view
        url = reverse("main:home")
        response = client.get(url,follow=True)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "main/home.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "home")

    def test_addGenome_view(self):
        # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)

        # get request to the validate view
        url = reverse("main:addGenome")
        response = client.get(url,follow=True)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "main/addGenome/addGenome.html")

        # Check the value of active_tab
        self.assertNotIn("active_tab", response.context)

    def test_genomeAdmin_view(self):
        # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)

        # get request to the validate view
        url = reverse("main:genomeAdmin")
        response = client.get(url,follow=True)
        # Check 200 status code
        self.assertEqual(response.status_code, 200)
        # Check template
        self.assertTemplateUsed(response, "main/admin/admin.html")
        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "admin")
        self.assertEqual(response.context["active_tab_admin"], "genome")
          
    def test_sequenceAdmin_view(self):
        # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)

        # get request to the validate view
        url = reverse("main:sequenceAdmin")
        response = client.get(url,follow=True)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        response = client.get(url)

        # Check template
        self.assertTemplateUsed(response, "main/admin/admin.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "admin")
        self.assertEqual(response.context["active_tab_admin"], "sequence")

    def test_accountAdmin_view(self):
        # Create a client
        admin_user = CustomUser.objects.create_superuser(email="super@user.com", firstName="admin",lastName="admin",role=3,password="ahhhhh")
        # Initialiser le client avec le super utilisateur
        client = Client()
        client.force_login(admin_user)

        # get request to the validate view
        url = reverse("main:accountAdmin")
        response = client.get(url,follow=True)

        # Check 200 status code
        self.assertEqual(response.status_code, 200)

        # Check template
        self.assertTemplateUsed(response, "main/admin/admin.html")

        # Check the value of active_tab
        self.assertEqual(response.context["active_tab"], "admin")
        self.assertEqual(response.context["active_tab_admin"], "account")
    


from .views import kind_of_sequence, blast
from unittest.mock import patch

class BlastTest(TestCase):
    def test_kind_of_sequence(self):
        self.assertEqual(kind_of_sequence("ATGC"), "nuc")
        self.assertEqual(kind_of_sequence("ATGCRYKMSWBDHVN"), "nuc")
        self.assertEqual(kind_of_sequence("ATcGaTt"), "nuc")
        self.assertEqual(kind_of_sequence("ACDEFGHIKLMNPQRSTVWY"), "prot")
        self.assertEqual(kind_of_sequence("acdefghiklmnpqrstvwy"), "prot")
        self.assertEqual(kind_of_sequence("invalid_sequence"), "pb_seq")
        self.assertEqual(kind_of_sequence("pbseq99"), "pb_seq")

    @patch('main.views.NCBIWWW.qblast')
    @patch('main.views.SearchIO.read')
    def test_blast_view(self, mock_read, mock_qblast):
        # Mocking the NCBIWWW.qblast and SearchIO.read functions
        mock_qblast.return_value = "fake_result_handle"
        mock_read.return_value = "fake_blast_results"

        url = reverse('main:blast')
        response = self.client.post(url, {'sequence': 'ATGC', 'program': 'blastn', 'alignments': 5})

        self.assertEqual(response.status_code, 200)
        self.assertEqual(mock_qblast.call_count, 1)
        self.assertEqual(mock_read.call_count, 1)
        self.assertIn('results', response.context)
        self.assertEqual(response.context['results'], 'fake_blast_results')

    def test_blast_view_error(self):
        url = reverse('main:blast')
        response = self.client.post(url, {'sequence': 'invalid_sequence', 'program': 'blastn', 'alignments': 5},follow=True)

        self.assertEqual(response.status_code, 200)
        self.assertIn('error_message', response.context)
        self.assertEqual(response.context['error_message'], 'Please verify that your query is a protein or a nuc sequence')
